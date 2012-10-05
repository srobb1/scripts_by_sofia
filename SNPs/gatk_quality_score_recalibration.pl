#!/usr/bin/perl -w
use strict;
use File::Spec;
##use realigned.bam reference.fa
my $ref = shift;
my @realigned_bams = @ARGV;

foreach my $bam (@realigned_bams) {
  my $bampath = File::Spec->rel2abs($bam);
  my @path    = split '/', $bampath;
  my $file    = pop @path;
  my ($root) = $file =~ /(.+)\.bam/;
  my $path = join '/', @path;

  my $name    = "$path/$root";
  my $cwd     = `pwd`;
  my $tempDir = '/scratch';
  open OUTSH, ">$name.gatk_QSR.sh"
    or die "Can't open $name.gatk_QSR.sh for writing\n";
  print OUTSH "#PBS -N $root\n";
  print OUTSH "hostname\n";
  print OUTSH "cd $cwd\n";
  print OUTSH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
  print OUTSH "
#if [ ! -f $name.mateFixed.bam ]; then
#/usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir  -jar picard/FixMateInformation.jar INPUT=$file OUTPUT=$name.mateFixed.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true
#fi

if [ ! -f $name.cc.csv ]; then
/usr/local/java/oracle-jdk1.7.0_01/bin/java -Djava.io.tmpdir=\$tmp_dir -Xmx6g -jar /srv/zpools/tern.ib_bigdata/home/stajichlab/shared/pkg/GATK/GenomeAnalysisTK.jar  -R $ref -I $file -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile $name.cc.csv
fi

if [ ! -f $name.recal.bam ]; then
/usr/local/java/oracle-jdk1.7.0_01/bin/java -Xmx6g -jar /srv/zpools/tern.ib_bigdata/home/stajichlab/shared/pkg/GATK/GenomeAnalysisTK.jar -R $ref -I $file -T TableRecalibration --out $name.recal.bam -recalFile $name.cc.csv
fi

rm -rf \$tmp_dir

";
}
