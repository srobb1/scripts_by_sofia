#!/usr/bin/perl -w
use strict;
use File::Spec;
my $GENOME = shift;
my $BASE = shift;
my @bams = @ARGV;
my $input;
#my $java = '/usr/local/java/oracle-jdk1.7.0_01/bin/java';
#my $GATK='/shared/stajichlab/pkg/GATK/GenomeAnalysisTK.jar';
my $tempDir = '/scratch';
$GENOME = File::Spec->rel2abs( $GENOME );
foreach my $file (@bams){ 
  my $path = File::Spec->rel2abs( $file );
  $input .=  "-I $path "
}
my $cwd = `pwd`;
open OUTSH ,">find_SNPs_VQSR.sh";
print OUTSH "
module load GATK
module load java/1.7.0_11

if [ -d $tempDir ]; then
 tmp_dir=`mktemp --tmpdir=$tempDir -d`
else
 tmp_dir=`mktemp --tmpdir=/tmp -d`
fi
cd $cwd

if [ ! -f  $BASE.raw.vcf ]; then
java -Djava.io.tmpdir=\$tmp_dir -Xmx128g -jar \$GATK -T UnifiedGenotyper -R $GENOME $input  -o $BASE.raw.vcf -glm SNP --read_filter BadCigar -nt 48 --metrics_file $BASE.info
fi


if [ ! -f  $BASE.VQSR.out ]; then
java -Xmx4g -jar \$GATK \\
   -T VariantRecalibrator \\
   -R $GENOME \\
   -input $BASE.raw.vcf \\
   -mode SNP \\
   -resource:dbsnp,known=true,training=true,truth=true,prior=6.0 $BASE.raw.vcf \\
   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ0 \\
   -recalFile $BASE.VQSR.out \\
   -tranchesFile $BASE.VQSR.tranches \\
   -rscriptFile $BASE.VQSR.plots.R \\
   --maxGaussians 4 \\
   --percentBadVariants 0.05
fi

if [ ! -f $BASE.VQSR.filtered.vcf ] ; then
java -Xmx3g -jar \$GATK \\
   -T ApplyRecalibration \\
   -R $GENOME \\
   -input $BASE.raw.vcf \\
   --ts_filter_level 99.0 \\
   -tranchesFile $BASE.VQSR.tranches \\
   -recalFile $BASE.VQSR.out \\
   -mode SNP \\
   -o $BASE.VQSR.filtered.vcf
fi



rm -rf \$tmp_dir
";


