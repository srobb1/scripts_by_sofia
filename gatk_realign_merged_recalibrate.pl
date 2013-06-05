#!/usr/bin/perl -w

##
## http://www.broadinstitute.org/gatk/guide/topic?name=best-practices
## Best: multi-sample realignment with known sites and recalibration
##
## for each sample
##    lanes.bam <- merged lane.bams for sample
##    dedup.bam <- MarkDuplicates(lanes.bam)
##    realigned.bam <- realign(dedup.bam) [with known sites included if available]
##    recal.bam <- recal(realigned.bam)
##    sample.bam <- recal.bam
##
##



use strict;
use File::Spec;
my $dir       = shift or die "dir genomefile strain specification\n";
my $genome    = shift;
my $knownSites = shift;
my $strain    = shift;
my $specification = shift;
my $change_readgroup = 1;
if (!defined $strain and !defined $specification){
  $change_readgroup = 0;
}


my $JAVA = '/usr/local/java/oracle-jdk1.7.0_01/bin/java';
#my $GATK='/shared/stajichlab/pkg/GATK/GenomeAnalysisTK.jar';
my $center = 'UCR';
my $type ='Genomic';
my $platform  = 'illumina';
$dir =~ s/\/$//;
my $path = File::Spec->rel2abs( $dir );
my @bam_files = <$dir/*bam>;
my $input;
foreach my $bam (@bam_files) {
  my $bampath = File::Spec->rel2abs( $bam );
  $input .= "INPUT=$bampath "; 
}
  my $name = "$strain.$specification";
  my $cwd = `pwd`;
  my $tempDir = '/scratch';
  open OUTSH , ">$name.gatk_realign.sh" or die "Can't open $name.gatk_realign.sh for writing\n";
  print OUTSH "#PBS -N $name\n";
  print OUTSH "#PBS -l nodes=1:ppn=8\n";
  print OUTSH "hostname\n";
  print OUTSH "cd $cwd\n";
  print OUTSH "
if [ -d $tempDir ]; then
 tmp_dir=`mktemp --tmpdir=$tempDir -d`
else
 tmp_dir=`mktemp --tmpdir=/tmp -d`
fi

module load GATK
module load picard

\n";
 
 print OUTSH "echo \"start merge\"\n";
  print OUTSH "
##
## lanes.bam <- merged lane.bams for sample
##

if [ ! -f $name.merged.bam ]; then
  $JAVA -Xmx24g -jar /usr/local/java/common/lib/picard-tools/MergeSamFiles.jar OUTPUT=$name.merged.bam SORT_ORDER=coordinate ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=true $input CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
fi
echo \"end merge\"
echo \"start sort\"
if [ ! -f $name.sort_ordered.bam ]; then
  $JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/ReorderSam.jar INPUT=$name.merged.bam OUTPUT=$name.sort_ordered.bam REFERENCE=$genome MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;
fi
echo \"end sort\"
";
if ($change_readgroup){
  print OUTSH "
echo \"start Readgroup\"
if [ ! -f $name.RG.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/AddOrReplaceReadGroups.jar INPUT=$name.sort_ordered.bam OUTPUT=$name.RG.bam SORT_ORDER=coordinate RGLB=$strain RGID=$strain RGSM=$specification RGPL=$platform RGCN=$center RGPU=$type VALIDATION_STRINGENCY=SILENT;
fi
echo \"end Readgroup\"
";
}else {
  print OUTSH "
  ln -s $name.sort_ordered.bam $name.RG.bam
"
}
print OUTSH "
##
## dedup.bam <- MarkDuplicates(lanes.bam)
##

echo \"start dedup\"
if [ ! -f $name.dedup.bam ]; then
  $JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/MarkDuplicates.jar INPUT=$name.RG.bam OUTPUT=$name.dedup.bam METRICS_FILE=$name.dedup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;
fi
echo \"end dedup\"

##
## realigned.bam <- realign(dedup.bam) [with known sites included if available]
##

echo \"start interval\"
if [ ! -f $name.intervals ]; then
  $JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar \$GATK -T RealignerTargetCreator -R $genome -o $name.intervals -I $name.dedup.bam --known $knownSites
fi
echo \"end interval\"


echo \"start realign\"
if [ ! -f $name.realign.bam ]; then
  $JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar \$GATK -T IndelRealigner -R $genome -targetIntervals $name.intervals -I $name.dedup.bam -o $name.realign.bam 
fi
echo \"end realign\"

##
##  recal.bam <- recal(realigned.bam)
##

echo \"start recalibate\"
if [ ! -f $name.recal_data.grp ] ; then 
 $JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx4g -jar \$GATK \\
 -T BaseRecalibrator \\
 -I $name.realign.bam \\
 -R $genome \\
 -knownSites $knownSites \\
 -o $name.recal_data.grp 
fi

if [ ! -f $name.recal.bam ] ; then 
$JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx4g -jar \$GATK \\
   -T PrintReads \\
   -R $genome \\
   -I $name.realign.bam \\
   -BQSR $name.recal_data.grp \\
   -o $name.recal.bam
fi
echo \"end recalibate\"

echo \"remove tmp dir\"

rm -rf \$tmp_dir

echo \"done\"
";

