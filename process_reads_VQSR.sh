#PBS -l nodes=1;ppn=8
GENOME=$1
FQ_1=$2
FQ_2=$3
BASE=$4
knownSites=$5
CWD=$6

#CWD=$PBS_O_WORKDIR
#if [ -z $CWD ] ; then
#  CWD=`pwd`
#fi

if [ $# -eq 0 ] # same effect as:  if [ -z "$1" ]
# $1 can exist, but be empty:  zmore "" arg2 arg3
then
  echo "Usage: `basename $0` GenomeFilePath FQ_1 FQ_2 SampleShortName VCFofKnownSites WorkingDir" >&2
  exit 1 
fi 


cd $CWD

# include some programs in our path on the biocluster system (would be different on other systems)
#module load bwa/0.6.2
module load samtools/0.1.18-r580
module load GATK
module load picard
module load sickle
module load bwa/0.6.2
module load java/1.7.0_11
# index the genome for alignment
if [ ! -f $GENOME.bwt ]; then
 bwa index -a bwtsw $GENOME
fi


if [ ! -d clean_fq ]; then
  mkdir $CWD/clean_fq
fi

CLEAN_1=$CWD/clean_fq/"$BASE"_p1.fq
CLEAN_2=$CWD/clean_fq/"$BASE"_p2.fq
CLEAN_U=$CWD/clean_fq/"$BASE"_unpaired.fq

# trim some reads before processing
# This strain is W303 of yeast, you can either use this short name or the original SRR567756 if you like
if [ ! -f $CLEAN_1 ]; then
 sickle pe -f $FQ_1 -r $FQ_2 -o $CLEAN_1 -p $CLEAN_2 -s $CLEAN_U -t sanger -q 20 -l 50
 echo "sickle pe -f $FQ_1 -r $FQ_2 -o $CLEAN_1 -p $CLEAN_2 -s $CLEAN_U -t sanger -q 20 -l 50"
fi

if [ ! -f $BASE.sam ]; then
 echo "bwa aln -t 8 -q 20 $GENOME $CLEAN_1 > ${BASE}_p1.sai" 
 bwa aln -t 8 -q 20 $GENOME $CLEAN_1 > "$BASE"_p1.sai 
 echo "bwa aln -t 8 -q 20 $GENOME $CLEAN_2 > "$BASE"_p2.sai"
 bwa aln -t 8 -q 20 $GENOME $CLEAN_2 > "$BASE"_p2.sai
 echo "bwa sampe $GENOME "$BASE"_p1.sai "$BASE"_p2.sai  $CLEAN_1 $CLEAN_2 | samtools view -Sb - -o $BASE.bam" 
 bwa sampe $GENOME "$BASE"_p1.sai "$BASE"_p2.sai  $CLEAN_1 $CLEAN_2 > $BASE.sam #| samtools view -Sb - -o $BASE.bam 
 #samtools index $BASE.bam
fi

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

#make the SAM file, then the BAM file as a sorted file
if [ ! -f $BASE.bam ]; then
 # now sort: ask for 3gb of memory in case this is big datafile
java -Xmx3g -jar $PICARD/SortSam.jar I=$BASE.sam O=$BASE.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE
 samtools flagstat $BASE.bam > $BASE.flagstat
# java -Xmx3g -jar $PICARD/ReorderSam.jar INPUT=$BASE.bam OUTPUT=$BASE.sort_ordered.bam SORT_ORDER=coordinate REFERENCE=$GENOME MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
fi

# Mark duplicate reads (usually where the forward and reverse are identical, indicating a
# PCR bias
if [ ! -f $BASE.dedup.bam ]; then
java -Xmx2g -jar $PICARD/MarkDuplicates.jar I=$BASE.bam \
  O=$BASE.dedup.bam METRICS_FILE=$BASE.dedup.metrics \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
 samtools flagstat $BASE.dedup.bam > $BASE.dedup.flagstat
fi

# Fix the ReadGroups - required by GATK
# right now the read groups aren't set in the depdup.bam file
if [ ! -f $BASE.RG.bam ]; then
 java -Xmx2g -jar $PICARD/AddOrReplaceReadGroups.jar I=$BASE.dedup.bam O=$BASE.RG.bam \
  SORT_ORDER=coordinate CREATE_INDEX=TRUE \
   RGID=$BASE RGLB=$BASE RGPL=Illumina RGPU=Genomic RGSM=$BASE \
   VALIDATION_STRINGENCY=SILENT
fi

# Identify where the variants are to realign around these
# this includes Indels
if [ ! -f $BASE.intervals ]; then
 java -Xmx2g -jar $GATK -T RealignerTargetCreator \
 -R $GENOME \
 -o $BASE.intervals \
 -I $BASE.RG.bam \
 --known $knownSites
fi

# realign the BAM file based on the intervals where there are polymorphism
if [ ! -f $BASE.realign.bam ]; then
 java -Xmx2g -jar $GATK -T IndelRealigner \
  -R $GENOME \
  -targetIntervals $BASE.intervals -I $BASE.RG.bam -o $BASE.realign.bam
fi


##
##  recal.bam <- recal(realigned.bam)
##

if [ ! -f $BASE.recal_data.grp ] ; then 
 java -Xmx6g -jar $GATK \
 -T BaseRecalibrator \
 -I $BASE.realign.bam \
 -R $GENOME \
 -knownSites $knownSites \
 -o $BASE.recal_data.grp \
 --plot_pdf_file $BASE.recal_data.pdf 
fi

if [ ! -f $BASE.recal.bam ] ; then 
java -Xmx4g -jar $GATK \
   -T PrintReads \
   -R $GENOME \
   -I $BASE.realign.bam \
   -BQSR $BASE.recal_data.grp \
   -o $BASE.recal.bam
fi


# Call the SNPs from this BAM file generating a VCF file
# using 4 threads (-nt 4) and only calling SNPs, INDELs could be call too
# with the -glm BOTH or -glm INDEL

## add call specific bases
if [ ! -f  $BASE.raw.vcf ]; then
java -Xmx3g -jar $GATK -T UnifiedGenotyper \
  -glm SNP \
  -I $BASE.recal.bam \
  -R $GENOME \
  -o $BASE.raw.vcf \
  -nt 4 \
  --dbsnp $knownSites \
  --genotyping_mode --genotyping_mode
fi

if [ ! -f  $BASE.VQSR.out ]; then
java -Xmx4g -jar $GATK \
   -T VariantRecalibrator \
   -R $GENOME \
   -input $BASE.raw.vcf \
   -mode SNP \
   -resource:dbsnp,known=true,training=true,truth=true,prior=6.0 $BASE.raw.vcf \
   -an QD -an HaplotypeScore -an MQRankSum -an ReadPosRankSum -an FS -an MQ0 \
   -recalFile $BASE.VQSR.out \
   -tranchesFile $BASE.VQSR.tranches \
   -rscriptFile $BASE.VQSR.plots.R \
   --maxGaussians 4 \
   --percentBadVariants 0.05
fi

if [ ! -f $BASE.VQSR.filtered.vcf ] ; then 
 java -Xmx3g -jar $GATK \
   -T ApplyRecalibration \
   -R $GENOME \
   -input $BASE.raw.vcf \
   --ts_filter_level 99.0 \
   -tranchesFile $BASE.VQSR.tranches \
   -recalFile $BASE.VQSR.out \
   -mode SNP \
   -o $BASE.VQSR.filtered.vcf
fi

module load vcftools

# run VCF tools to convert the filtered VCF file into tab-delimited
# for some simple look at the SNPs
# would also do other work with the VCF file in vcftools to look at summary statistics
vcf-to-tab < $BASE.VQSR.filtered.vcf > $BASE.VQSR.filtered.tab

