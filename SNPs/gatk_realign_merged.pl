#!/usr/bin/perl -w
use strict;
use File::Spec;
my $dir       = shift or die "dir genomefile strain specification\n";
my $genome    = shift;
my $strain    = shift;
my $specification = shift;
my $change_readgroup = 1;
if (!defined $strain and !defined $specification){
  $change_readgroup = 0;
}
my $GATK='/shared/stajichlab/pkg/GATK/GenomeAnalysisTK.jar';
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
  print OUTSH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
  print OUTSH "echo \"start merge\"\n";
  print OUTSH "
if [ ! -f $name.merged.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Xmx24g -jar /usr/local/java/common/lib/picard-tools/MergeSamFiles.jar OUTPUT=$name.merged.bam SORT_ORDER=coordinate ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=true $input CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
fi
echo \"end merge\"
echo \"start sort\"
if [ ! -f $name.sort_ordered.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/ReorderSam.jar INPUT=$name.merged.bam OUTPUT=$name.sort_ordered.bam REFERENCE=$genome MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;
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
echo \"start dedup\"
if [ ! -f $name.dedup.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/MarkDuplicates.jar INPUT=$name.RG.bam OUTPUT=$name.dedup.bam METRICS_FILE=$name.dedup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;
fi
echo \"end dedup\"
echo \"start interval\"
if [ ! -f $name.intervals ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar $GATK -T RealignerTargetCreator -R $genome -o $name.intervals -I $name.dedup.bam;
fi
echo \"end interval\"
echo \"start realign\"
if [ ! -f $name.realign.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar $GATK -T IndelRealigner -R $genome -targetIntervals $name.intervals -I $name.dedup.bam -o $name.realign.bam;
fi
echo \"end realign\"
echo \"remove tmp dir\"

rm -rf \$tmp_dir
echo \"done\"
";


#  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/SortSam.jar INPUT=$bampath OUTPUT=$name.sort_ordered.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT;


