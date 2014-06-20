cd /rhome/robb/project/bisulfide/mPing_primers

#FILE=$1
FILE=A119_2.mping.confident_nonref_genomeflanks.fa
BASE=`basename $FILE .confident_nonref_genomeflanks.fa`
FASTA=$BASE.fa
perl -p -e 's/\>(\S+) (\S+) (\S+)/\>$3 mPing\@$1/' $FILE > $FASTA
perl ~/project/bisulfide/mPing_primers/makePrimers.pl -f $FASTA -optimalProductLength 150 -seqTarget 98,3 -primerLenOpt 30 -primerLenMin 27 -primerLenMax 34 -tmMin 62 -tmOpt 65 -tmMax 70 > $BASE.primers.txt
perl get_methylation_on_primers_faster.pl $BASE.primers.txt NB   > $BASE.primers.NBmethylation.txt
perl get_methylation_on_primers_faster.pl $BASE.primers.txt A119 > $BASE.primers.A119methylation.txt
#perl get_methylation_on_primers.pl $BASE.primers.txt NB   > $BASE.primers.NBmethylation.txt
#perl get_methylation_on_primers.pl $BASE.primers.txt A119 > $BASE.primers.A119methylation.txt
