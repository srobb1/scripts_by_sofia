#!/usr/bin/perl -w
use strict;
use File::Spec;
my $dir       = shift;
my $genome    = shift;
my $strain    = shift;
my $specification = shift;

if (!defined $dir or !defined $genome or !defined $strain or !defined $specification){
  die "gatk_realign.pl dir_of_bamfiles path_2_genome strain strain_specificatio
ex gatk_realign.pl ~/rice/SNPs/A123_0 ~/Wessler-Rice/Genome/index/MSUr7.all.fa A123 A123_0
"
}

my $center = 'UCR';
my $type ='Genomic';
my $platform  = 'illumina';
$dir =~ s/\/$//;
my @bam_files = <$dir/*bam>;
foreach my $bam (@bam_files) {
  my $bampath = File::Spec->rel2abs( $bam );
  my @path = split '/', $bampath;
  my $file = pop @path;
  my ($root) = $file =~ /(.+)\.bam/;
  my ($lib) = $root =~ /(FC\d+_\d+)/;
  my $path = join '/' , @path;

  my $name = "$path/$root";
  my $cwd = `pwd`;
  my $tempDir = '/scratch';
  open OUTSH , ">$name.gatk_realign.sh" or die "Can't open $name.gatk_realign.sh for writing\n";
  print OUTSH "#PBS -N $root\n";
  print OUTSH "hostname\n";
  print OUTSH "cd $cwd\n";
  print OUTSH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
  print OUTSH "
if [ ! -f $name.sort_ordered.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/ReorderSam.jar INPUT=$bampath OUTPUT=$name.sort_ordered.bam REFERENCE=$genome MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;
fi
if [ ! -f $name.RG.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/AddOrReplaceReadGroups.jar INPUT=$name.sort_ordered.bam OUTPUT=$name.RG.bam SORT_ORDER=coordinate  RGLB=$strain RGID=$strain RGSM=$specification RGPL=$platform RGCN=$center RGPU=$type VALIDATION_STRINGENCY=SILENT;
fi
if [ ! -f $name.dedup.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/MarkDuplicates.jar INPUT=$name.RG.bam OUTPUT=$name.dedup.bam METRICS_FILE=$name.dedup.metrics MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT;
fi
if [ ! -f $name.intervals ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /srv/zpools/tern.ib_bigdata/home/stajichlab/shared/pkg/GATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $genome -o $name.intervals -I $name.dedup.bam;
fi
if [ ! -f $name.realign.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /srv/zpools/tern.ib_bigdata/home/stajichlab/shared/pkg/GATK/GenomeAnalysisTK.jar -T IndelRealigner -R $genome -targetIntervals $name.intervals -I $name.dedup.bam -o $name.realign.bam;
fi

rm -rf \$tmp_dir

";
}

#  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/SortSam.jar INPUT=$bampath OUTPUT=$name.sort_ordered.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT;


#/usr/local/java/oracle-jdk1.7.0_01/bin/java -Xmx24g -jar /usr/local/java/common/lib/picard-tools/MergeSamFiles.jar OUTPUT=/dev/shm/HEG4.all.bam SORT_ORDER=coordinate ASSUME_SORTED=true MERGE_SEQUENCE_DICTIONARIES=true USE_THREADING=true INPUT=3kb/rice3kb_FC70_5.realign.HEG4.bam INPUT=500bp/rice500bp_FC52_7.realign.HEG4.bam INPUT=500bp/rice500bp_FC52_8.realign.HEG4.bam INPUT=500bp/rice500bp_FC67_1.realign.HEG4.bam INPUT=500bp/rice500bp_FC67_2.realign.HEG4.bam INPUT=500bp/rice500bp_FC70_1.realign.HEG4.bam INPUT=500bp/rice500bp_FC70_2.realign.HEG4.bam INPUT=500bp/rice500bp_FC70_3.realign.HEG4.bam INPUT=5kb/rice5kb_FC20_3.realign.HEG4.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
