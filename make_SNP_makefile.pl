#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use File::Spec;
use FindBin qw($RealBin);

my $SCRIPTS     = $RealBin;
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

if ( !defined @ARGV ) {
  &getHelp();
}
my $DESTINATION_DIR;
my $FASTQ_DIR;
my $RENAME_FASTQ_DIR;
my $BARCODE_SAMPLE_FILE;
my $QUEUE_Q;
my $QUEUE_L;

my $dbSNP  = '/shared/wesslerlab/Rice/HEG4/dbSNP/dbSNPs_VQSR.vcf';
my $GENOME = '/shared/wesslerlab/Rice/Genome/index/MSU_r7.corrected.fa';

GetOptions(
  'o|output_dir:s'     => \$DESTINATION_DIR,
  'f|fastq_dir:s'      => \$FASTQ_DIR,
  'r|rename_fastq_dir:s'  => \$RENAME_FASTQ_DIR,
  'i|barcode_sample:s' => \$BARCODE_SAMPLE_FILE,
  'd|dbSNPs:s'         => \$dbSNP,
  'g|genome:s'         => \$GENOME,
  'q|qsub_q:s'         => \$QUEUE_Q,
  'l|qsub_1:s'         => \$QUEUE_L,
  'h|help'             => \&getHelp,
);

sub getHelp {
  print " 
Assistants in creating the config file for makefile to process reads for SNP calling.

usage: $0 -o output_dir -f fastq_dir -i barcode_sample_id_file -q qsub_options -h this message

  o|output_dir:s      	=> name of directory for all outputed directories and files
  f|fastq_dir:s       	=> name of directory of fastq files to be processed
  r|rename_fastq_dir:s 	=> name of directory of in which the renamed fastq files will be placed
  i|barcode_sample:s  	=> name of file with Sample_names and illumina barcode ids
  g|genome:s            => name of genome file 
  d|dbSNPs:s            => name of dbSNP file 
  q|qsub_q:s            => qsub -q options, if a queue is available. ex: -q highmem
  l|qsub_l:s            => qsub -l options, if a queue is available. ex: -l nodes=1:ppn=8
  h|help 		=> this message

";
  exit 1;
}

my $outdir;
if ( defined $DESTINATION_DIR and -e $DESTINATION_DIR ) {
  $outdir = File::Spec->rel2abs($DESTINATION_DIR);
}
else {
  print
"Please provide the directory for output. $DESTINATION_DIR is not valid or does not exist\n";
  getHelp();
}
my $fastq_dir;
if ( defined $FASTQ_DIR and -e $FASTQ_DIR ) {
  $fastq_dir = File::Spec->rel2abs($FASTQ_DIR);
}
else {
  print
"Please provide the directory in which your FASTQ files can be found. $FASTQ_DIR is not valid\n";
  getHelp();
}
my $rename_fastq_dir;
if ( defined $RENAME_FASTQ_DIR and -e $RENAME_FASTQ_DIR ) {
  $rename_fastq_dir = File::Spec->rel2abs($RENAME_FASTQ_DIR);
}
else {
  print
"Please provide the directory in which your FASTQ files can be found. $RENAME_FASTQ_DIR is not valid\n";
  getHelp();
}
my $id_file;
if ( defined $BARCODE_SAMPLE_FILE and -e $BARCODE_SAMPLE_FILE ) {
  $id_file = File::Spec->rel2abs($BARCODE_SAMPLE_FILE);
}
else {
  print
"Please provide the path to the file of Sample names and Illumina barcode IDs. $BARCODE_SAMPLE_FILE is not valid or does not exist\n";
  getHelp();
}
my $dbSNP_file;
if ( defined $dbSNP and -e $dbSNP ) {
  $dbSNP_file = File::Spec->rel2abs($dbSNP);
}
else {
  print
"Please provide the path to the file of Genotype Postitions and alleles. $dbSNP is not valid or does not exist\n";
  getHelp();
}
my $genome_file;
my $genome_base;
if ( defined $GENOME and -e $GENOME ) {
  $genome_file = File::Spec->rel2abs($GENOME);
  my @genome = split /\//, $genome_file;
  my $file = pop @genome;
  ($genome_base) = $file =~ /(.+)\.(fa|fasta)$/;
}
else {
  print
"Please provide the path to the Genome file. $GENOME is not valid or does not exist\n";
  getHelp();
}

my $qsub_q = defined $QUEUE_Q ? "-q $QUEUE_Q" : '';
my $qsub_l = defined $QUEUE_L ? "-l $QUEUE_L" : '';

open MAKEFILE, ">$fastq_dir/Makefile"
  or die "Can't open $fastq_dir/Makefile for writing $!\n";

print MAKEFILE "GENOME=$genome_file
HOME=$current_dir
DESTINATION=$outdir
RENAME_DIR=$rename_fastq_dir
OUTDIR=\$(DESTINATION)/$genome_base
SCRIPTS=$SCRIPTS
FINAL_RENAME=\$(HOME)/final_SNP_renamed.txt
dbSNP=$dbSNP_file

all: SNP.todo.txt
clean: 
	rm \$(FINAL_RENAME)
	rm \$(HOME)/SNPcalling.todo.txt 
	
SNPcalling.todo.txt: \$(FINAL_RENAME)
	for i in `cat \$^` ; do\\
	 SHORT=`echo \$\$i|cut -d ',' -f1`;\\
	 FQ_1=`echo \$\$i|cut -d ',' -f2`;\\
	 FQ_2=`echo \$\$i|cut -d ',' -f3`;\\
	 printf \"\\
	 date\\n\\
	 if [ ! -d \$(OUTDIR) ]; then\\n\\
	  mkdir -p \$(OUTDIR)\\n\\
	 fi\\n\\
	 \$(SCRIPTS)/process_reads_SNPs.sh \$(GENOME) \$\$FQ_1 \$\$FQ_2 \$\$SHORT \$(dbSNP) \$(OUTDIR)\\n\\
	 date\\n\" > \$(HOME)/run.SNPcalling.\$\$SHORT.sh;\\
	 JOB=`qsub $qsub_q $qsub_l \$(HOME)/run.SNPcalling.\$\$SHORT.sh`;\\
	 echo \"`date` \$\$JOB run.SNPcalling.\$\$SHORT.sh\" >> \$\@;\\
	done;\\
	printf \"\\nThe following jobs have been submitted to the queue:\\n;\"\\
	cat  \$\@;

\$(FINAL_RENAME): \$(SCRIPTS)/rename_files.pl \$(RENAME_DIR) \$(HOME) $id_file
	\$^ > \$\@

.PHONY: all clean

";

#`$SCRIPTS/Makefile`;
print "Now run $SCRIPTS/Makefile like this:
cd $SCRIPTS
make all
";
## this script will go thru each entry in the sample/illumina_barcode_id file (see below for file format)
## and rename the flowcell193_lane5_pair1_ATCACG.fastq formatted files to this format: RIL13_0_TTAGGC_FC193L5_p2.fq
## it will then create shell scripts which get submitted to the queue for trimming and genotyping

############## input file formats ##############################################
####### $BARCODE_SAMPLE_FILE format: #######
##cat FC193_RIL_12-16_43/name_illumina_id.txt
##RIL12_0 27
##RIL13_0 3
##RIL16_0 22

####### final_trim.txt format: ##############
##cat final_trim.txt
##RIL13_0_TTAGGC_FC193L5,/rhome/robb/Wessler-Rice/RIL/Illumina/RIL13_0/RIL13_0_TTAGGC_FC193L5_p1.fq,/rhome/robb/Wessler-Rice/RIL/Illumina/RIL13_0/RIL13_0_TTAGGC_FC193L5_p2.fq
##RIL45_0_ATCACG_FC193L5,/rhome/robb/Wessler-Rice/RIL/Illumina/RIL45_0/RIL45_0_ATCACG_FC193L5_p1.fq,/rhome/robb/Wessler-Rice/RIL/Illumina/RIL45_0/RIL45_0_ATCACG_FC193L5_p2.fq

################################################################################

