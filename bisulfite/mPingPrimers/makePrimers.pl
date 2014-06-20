#!/usr/bin/perl -w

# this script will write a primer3 input file and will run parsePrimer3.pl 
# runPrimer3.pl -f fasta.nt  (Note -f is required)
# input a fasta file
# output a tab delimited file will all possible primers

use strict;
use Getopt::Long;
use Bio::SeqIO;
use IO::String;


my $seqTarget = '' ;
my ($gcOpt, $gcMin, $gcMax) =  (50,45,60) ;
my ($tmOpt, $tmMin, $tmMax)  = (55,50,65) ;
my ($primerLenOpt, $primerLenMin, $primerLenMax) = (20,18,25) ;
my $gcClamp = 1;
my $optimalProductLength = 'long';
my $file ='';
my $verbose = 0;
sub init(){
GetOptions (
              "seqTarget=s" => \$seqTarget,
              "gcOpt=i"           => \$gcOpt,
              "gcMin=i"           => \$gcMin,
              "gcMax=i"           => \$gcMax,
              "tmOpt=i"           => \$tmOpt,
              "tmMin=i"           => \$tmMin,
              "tmMax=i"           => \$tmMax,
              "primerLenOpt=i"           => \$primerLenOpt,
              "primerLenMin=i"           => \$primerLenMin,
              "primerLenMax=i"           => \$primerLenMax,
              "optimalProductLength=i"   => \$optimalProductLength,
              "gcClamp=i"         => \$gcClamp,
              "f|file=s"            => \$file,      # string
              "verbose"  => \$verbose);   # flag

warn qq(
   "seqTarget=s" => $seqTarget,
              "gcOpt=i"           => $gcOpt,
              "gcMin=i"           => $gcMin,
              "gcMax=i"           => $gcMax,
              "tmOpt=i"           => $tmOpt,
              "tmMin=i"           => $tmMin,
              "tmMax=i"           => $tmMax,
              "primerLenOpt=i"    => $primerLenOpt,
              "primerLenMin=i"    => $primerLenMin,
              "primerLenMax=i"    => $primerLenMax,
              "optimalProductLength=i"   => $optimalProductLength,
              "gcClamp=i"         => $gcClamp,
              "f|file=s"          => $file,  
              "verbose"           => $verbose 
);

}

sub usage() {
        print STDERR <<EOF;

To run this program  ...

usage: $0 [-h] [-g GC-opt,min,max] [-t Tm-opt,min,max][-l PrimerLen-opt,min,max][-f file]

required:
        -f file : fasta file
optional:
        -h      : this (help) message
        -g      : GC content -- optimal,min,max (default is 50,45,60)
        -t      : Melting Temperature -- optimal,min,max (default is 55,50,65)
        -l      : Primer Length  -- optimal min max (default is 20,18,25)

        example: $0 -f seqs.fasta
        example: $0 -g 50,45,60 -t 60,55,65 -f seqs.fasta


EOF
        exit;
}

init();

####################
#
# tests
#
####################


#test to see if file exists
die "\n\n***Oops the file: \'$file\' does not exist.\n Is it spelled right??\n\n" unless -e $file;

#replace carriage returns for newlines
`perl -p -i -e 's/\r/\n/g' $file`;

#test for fasta format
my $found = `grep ">" $file`;

if (!$found){
        die "\n\n***ERROR***
                Input file not in FASTA format.
                No primers were produced.

        EXAMPLE FASTA FILE

        >SEQNAME
        GTATATATGGGTAAATAAATTGAAAATTACTTGTCTTAT
        TCCTTCAGTAATTGCCGTAACTGAAAATTTTGATGCTGA
        ATAGAAATGAATTGAAGAAAAGCTGAAAACCTTGTGACC


"
;
}


####################
#
#  1.parsing of the required fasta file
#  2.writing the primer3 input file
#
####################
my $seqIO = Bio::SeqIO-> new(-file     => $file,
                             -format => 'fasta');
open OUTFILE,  ">$file.in_primer3" or die "Can't open $file.in_primer3";
my $count = 0;
while (my $seq_obj = $seqIO->next_seq) {
        my $id  = $seq_obj->id;
        my $desc = $seq_obj->desc;
        if (defined $desc){
          $id .= "|$desc";
        } 
        my $seq = $seq_obj->seq;
        my $total_len = $seq_obj->length;
        my $overlap='';
        if ($seq =~ /\-/){
          $total_len = 0;
          my @overlap;
          my @frags = split '-' , $seq;
          foreach my $frag (@frags){
            $total_len += length $frag;
            push @overlap, $total_len;
          }
          pop @overlap;
          $overlap = "\nSEQUENCE_OVERLAP_JUNCTION_LIST=".(join ' ' , @overlap);
          $seq = join ('',@frags);
        }
        my $len = length($seq);
	if ($optimalProductLength !~ /long/){
		$len = $optimalProductLength;

	}
        #my $len = length($seq);
        my $len_75 = int ($len * .75);
        my $len_50 = int ($len * .50);
        my $len_25 = int ($len * .25);
        if ($primerLenMax > $total_len){
		warn "** WARNING ** $id sequence is too short.  No primers will be made for this sequence.\n";
		next if $len_25 <= $primerLenMax;
	}	
	

        print OUTFILE

"PRIMER_SEQUENCE_ID=$id
TARGET=$seqTarget
SEQUENCE=$seq$overlap
PRIMER_PRODUCT_SIZE_RANGE=$len_75-$len $len_50-$len $len_25-$len 
PRIMER_GC_CLAMP=$gcClamp       
PRIMER_MIN_GC=$gcMin
PRIMER_OPT_GC_PERCENT=$gcOpt
PRIMER_MAX_GC=$gcMax
PRIMER_OPT_SIZE=$primerLenOpt
PRIMER_MIN_SIZE=$primerLenMin
PRIMER_MAX_SIZE=$primerLenMax
PRIMER_OPT_TM=$tmOpt
PRIMER_MIN_TM=$tmMin
PRIMER_MAX_TM=$tmMax
=\n";
	$count++
}



if ($count){ #proceed only if 1 or more appropriate sequences are provided

####################
#
# running the primer3
#
###################

	warn "Primers for $count seqeunce(s).\n\n";

	`cat $file.in_primer3  | primer3_core > $file.out_primer3`;

	#print "Your Primers are in $file.primers.txt\n\n";



####################
#
# running the parser
#
###################

	my $output = `./parsePrimer3.pl $file.out_primer3` ;

	print $output;
}

