GENOME_input=$1
FQ1_input=$2
FQ2_input=$3
BASE_input=$4
PWD=`pwd`
echo "
#PBS -q highmem

module load bowtie2/2.1.0
module load bismark
module load sickle


#### FILL OUT Variables below ########
CWD=$PWD
OUTPUT_DIR=\$CWD/output
GENOME=$GENOME_input
FQ1=\$CWD/reads/$FQ1_input
FQ2=\$CWD/reads/$FQ2_input
BASE=$BASE_input
##### FINISHED with VARIABLES ########


BASE1=`basename $FQ1_input .fq`
BASE2=`basename $FQ2_input .fq`


cd \$CWD

if [ ! -e \$OUTPUT_DIR ] ; then
  mkdir \$OUTPUT_DIR
fi

if [ ! -d clean_fq ]; then
  mkdir \$CWD/clean_fq
fi

CLEAN_1=\$CWD/clean_fq/\${BASE1}.fq
CLEAN_2=\$CWD/clean_fq/\${BASE2}.fq
CLEAN_U=\$CWD/clean_fq/\${BASE}_unpaired.fq

if [ ! -f \$CLEAN_1 ]; then
 sickle pe -f \$FQ1 -r \$FQ2 -o \$CLEAN_1 -p \$CLEAN_2 -s \$CLEAN_U -t sanger -q 20 -l 50
 echo \"sickle pe -f \$FQ1 -r \$FQ2 -o \$CLEAN_1 -p \$CLEAN_2 -s \$CLEAN_U -t sanger -q 20 -l 50\"
fi



if [ ! -e "\$GENOME/Bisulfite_Genome/CT_conversion/BS_CT.1.bt2" ] ; then
  bismark_genome_preparation \$GENOME --bowtie2 --verbose 2&> bismark.info
fi 

##USAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}
TMP_DIR=\$OUTPUT_DIR/tmp
bismark \$GENOME -1 \$CLEAN_1 -2 \$CLEAN_2 \$CLEAN_U --bowtie2 --prefix \$BASE --output_dir \$OUTPUT_DIR --temp_dir \$TMP_DIR -N 1 
echo \"bismark \$GENOME -1 \$CLEAN_1 -2 \$CLEAN_2 \$CLEAN_U --bowtie2 --prefix \$BASE --output_dir \$OUTPUT_DIR --temp_dir \$TMP_DIR -N 1\"
"
