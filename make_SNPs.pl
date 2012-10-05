#!/usr/bin/perl -w
use strict;
use File::Spec;
my $genome = shift;
my $descr = shift;
my @bams = @ARGV;
my $input;
my $java = '/usr/local/java/oracle-jdk1.7.0_01/bin/java';
my $GATK='/shared/stajichlab/pkg/GATK/GenomeAnalysisTK.jar';
my $tempDir = '/scratch';
$genome = File::Spec->rel2abs( $genome );
foreach my $file (@bams){ 
  my $path = File::Spec->rel2abs( $file );
  $input .=  "-I $path "
}
my $cwd = `pwd`;
open OUTSH ,">find_SNPs.sh";
print OUTSH "
tmp_dir=`mktemp --tmpdir=$tempDir -d`
cd $cwd
$java  -Djava.io.tmpdir=\$tmp_dir -Xmx128g -jar $GATK -T UnifiedGenotyper -R $genome $input  -o $descr.raw.vcf -glm SNP --read_filter BadCigar -nt 48 --metrics_file $descr.info
rm -rf \$tmp_dir
";


