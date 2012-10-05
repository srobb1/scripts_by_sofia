#!/usr/bin/perl -w
use strict;
use  File::Spec;
my $dir = shift;
my $strain = shift;
my $specification = shift;
my $center = 'UCR';
my $type ='Genomic';
my $platform  = 'illumina';
$dir =~ s/\/$//;
my $path = File::Spec->rel2abs( $dir );
my @bam_files = <$dir/*bam>;
my $input;
my $cwd = `pwd`;
foreach my $bam (@bam_files) {
  my ($base) = $bam =~ /(.+)\.bam$/; 
  my $bampath = File::Spec->rel2abs( $bam );
  my $name = "$strain.$specification";
  my $tempDir = '/scratch';
  `mkdir -p $cwd/RG.$name`;
  #Can't open RIL39.RIL39_1.changeReadGroup.sh for writing
  open OUTSH , ">$cwd/RG.$name/$bam.changeRG.sh" or die "Can't open $name.changeRG.sh for writing\n";
  print OUTSH "#PBS -N $bam.changeRG\n";
  print OUTSH "#PBS -l nodes=1:ppn=8\n";
  print OUTSH "hostname\n";
  print OUTSH "cd $cwd\n";
  print OUTSH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
  print OUTSH "
if [ ! -f $name.RG.bam ]; then
  /usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /usr/local/java/common/lib/picard-tools/AddOrReplaceReadGroups.jar INPUT=$bam OUTPUT=$cwd/RG.$name/$base.RG.bam SORT_ORDER=coordinate RGLB=$strain RGID=$strain RGSM=$specification RGPL=$platform RGCN=$center RGPU=$type VALIDATION_STRINGENCY=SILENT;
fi
";
}
