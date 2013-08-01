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
module load GATK/2.4-3-g2a7af43
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
echo "Sorting $BASE.sam and creating $BASE.bam"
if [ ! -f $BASE.bam ]; then
 # now sort: ask for 3gb of memory in case this is big datafile
  java -Xmx3g -jar $PICARD/SortSam.jar I=$BASE.sam O=$BASE.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE
  echo "getting stats on $BASE.bam"
  samtools flagstat $BASE.bam > $BASE.flagstat
# java -Xmx3g -jar $PICARD/ReorderSam.jar INPUT=$BASE.bam OUTPUT=$BASE.sort_ordered.bam SORT_ORDER=coordinate REFERENCE=$GENOME MAX_RECORDS_IN_RAM=1000000 CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
fi

# Mark duplicate reads (usually where the forward and reverse are identical, indicating a
# PCR bias
echo "Marking duplicate reads: Generating $BASE.dedup.bam"
if [ ! -f $BASE.dedup.bam ]; then
java -Xmx2g -jar $PICARD/MarkDuplicates.jar I=$BASE.bam \
  O=$BASE.dedup.bam METRICS_FILE=$BASE.dedup.metrics \
  CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
 samtools flagstat $BASE.dedup.bam > $BASE.dedup.flagstat
fi

# Fix the ReadGroups - required by GATK
# right now the read groups aren't set in the depdup.bam file
echo "Generating $BASE.RG.bam"
if [ ! -f $BASE.RG.bam ]; then
 java -Xmx2g -jar $PICARD/AddOrReplaceReadGroups.jar I=$BASE.dedup.bam O=$BASE.RG.bam \
  SORT_ORDER=coordinate CREATE_INDEX=TRUE \
   RGID=$BASE RGLB=$BASE RGPL=Illumina RGPU=Genomic RGSM=$BASE \
   VALIDATION_STRINGENCY=SILENT
fi

echo "Generating $BASE.intervals from $BASE.RG.bam"
# Identify where the variants are to realign around these
# this includes Indels
if [ ! -f $BASE.intervals ]; then
 java -Xmx2g -jar $GATK -T RealignerTargetCreator \
 -R $GENOME \
 -o $BASE.intervals \
 -I $BASE.RG.bam \
 --known $knownSites
fi

echo "Generating $BASE.realign.bam from $BASE.RG.bam"
# realign the BAM file based on the intervals where there are polymorphism
if [ ! -f $BASE.realign.bam ]; then
 java -Xmx2g -jar $GATK -T IndelRealigner \
  -R $GENOME \
  -targetIntervals $BASE.intervals -I $BASE.RG.bam -o $BASE.realign.bam
fi


##
##  recal.bam <- recal(realigned.bam)
##

echo "Generating $BASE.recal_data.grp from $BASE.realign.bam"
if [ ! -f $BASE.recal_data.grp ] ; then 
 java -Xmx4g -jar $GATK \
 -T BaseRecalibrator \
 -I $BASE.realign.bam \
 -R $GENOME \
 -knownSites $knownSites \
 -o $BASE.recal_data.grp 
fi
echo "Generating $BASE.recal.bam from $BASE.realign.bam"
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
echo "Generating VCF from $BASE.recal.bam"
## add call specific bases
if [ ! -f  $BASE.genotype.vcf ]; then
java -Xmx3g -jar $GATK -T UnifiedGenotyper \
  -glm SNP \
  -I $BASE.recal.bam \
  -R $GENOME \
  -o $BASE.genotype.vcf \
  -nt 4 \
  --dbsnp $knownSites \
  --genotyping_mode GENOTYPE_GIVEN_ALLELES \
  --alleles $knownSites \
  --output_mode EMIT_ALL_SITES
fi


module load vcftools

# run VCF tools to convert the VCF file into tab-delimited
# for some simple look at the Genotypes
# would also do other work with the VCF file in vcftools to look at summary statistics
echo "Converting vcf to tab"
vcf-to-tab < $BASE.genotype.vcf > $BASE.genotype.tab


## clean up
FILESIZE=$(stat -c%s "$BASE.genotype.vcf")
if [[ $FILESIZE > 5000 ]]; then
  echo "Cleaning up: removing files"
  rm $BASE.realign.bai
  rm $BASE.realign.bam
  rm $BASE.dedup.metrics
  rm $BASE.dedup.bai
  rm $BASE.dedup.bam
  rm $BASE.recal_data.grp
  rm $BASE.bam
  rm $BASE.bai
  rm $BASE.sam
  rm $BASE.intervals
  rm ${BASE}_p2.sai
  rm ${BASE}_p1.sai
  rm $BASE.RG.bai
  rm $BASE.RG.bam
fi
