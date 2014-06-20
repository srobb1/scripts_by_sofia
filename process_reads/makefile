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





clean: 
	rm $FASTQ_DIR/final_renamed.txt

 

genotype: final_renamed.txt
        ARGS=( `cat $^` )
	for i in "${ARGS[@]}"
	do
	 echo "
 	 date

 	 SHORT=`echo $i|cut -d "," -f1`
	 FQ_1 =`echo $i|cut -d "," -f2`
	 FQ_2 =`echo $i|cut -d "," -f3`

  	 CWD=\$PBS_O_WORKDIR
  	 if [ -z \$CWD ] ; then
    	  CWD=\`pwd\`
	 fi
  	 cd \$CWD
         if [ ! -d \$CWD/MSU_r7.corrected ]; then
          mkdir -p \$CWD/MSU_r7.corrected
         fi
         cd \$CWD/MSU_r7.corrected
         $SCRIPTS/process_reads_genotype.sh \$GENOME \$FQ_1 \$FQ_2 \$SHORT \$dbSNP \$CWD/MSU_r7.corrected 
      
         date
         " > $CWD/run.genotype.${SHORT}.sh

         qsub $Q $CWD/run.genotype.${SHORT}.sh

	done


final_renamed.txt: $SCRIPTS/rename_files.pl $DESTINATION_DIR $FASTQ_DIR $BARCODE_SAMPLE_FILE
	cd $FASTQ_DIR
	$^ > $@


set_variables:  process_reads.config
	CONFIG=( `cat` process_reads.config )
 	SCRIPTS=${CONFIG[0]}
	CWD=${CONFIG[1]}
	DESTINATION_DIR=${CONFIG[2]}
	FASTQ_DIR=${CONFIG[3]}
	BARCODE_SAMPLE_FILE=${CONFIG[4]}
	dbSNP=${CONFIG[5]}
	GENOME=${CONFIG[6]}
	QUEUE=${CONFIG[7]}
