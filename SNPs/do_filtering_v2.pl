#!/usr/bin/perl -w
use strict;
##added --filterExpression \"DP<10\" --filterName LowCoverage
my $vcf = shift;
my $REF = shift; ##reference fasta

my ($base) = $vcf =~ /(.+).vcf/;

my $JAVA='/home_stajichlab/jstajich/bigdata-shared/pkg/jre1.6.0_29/bin/java'; 
my $GATK='/opt/stajichlab/GATK/latest/GenomeAnalysisTK.jar';

my $tmpDir = '/scratch';
open SH, ">$base.filter.sh";

print SH "
tmp_dir=`mktemp --tmpdir=$tmpDir -d`
cd `pwd`;

$JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx64g -jar $GATK -R $REF -T VariantFiltration -o $base.filter_v2.vcf --variant $vcf --clusterWindowSize 10  --filterExpression \"QD<5.0\" --filterName QualByDepth --filterExpression \"HRun>=4\" --filterName HomopolymerRun --filterExpression \"QUAL < 60\" --filterName QScore --filterExpression \"MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)\" --filterName MapQual --filterExpression \"FS > 60.0\" --filterName FisherStrandBias --filterExpression \"HaplotypeScore > 13.0\" --filterName HaplotypeScore --filterExpression \"MQRankSum < -12.5\" --filterName MQRankSum --filterExpression \"ReadPosRankSum < -8.0\" --filterName ReadPosRankSum --filterExpression \"DP<10\" --filterName LowCoverage 


$JAVA -Djava.io.tmpdir=\$tmp_dir -Xmx64g -jar $GATK -R $REF -T SelectVariants --variant $base.filter_v2.vcf -o $base.selectedSNPs_v2.vcf --excludeFiltered
rm -rf \$tmp_dir
";
