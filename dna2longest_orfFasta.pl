#!/usr/bin/perl -w

use Bio::SeqIO;
use strict;
use Data::Dumper;

######################################################################
#
#this first section will translate each fasta entry into 6 reading
#frames. The translated seq will be printed to a file called
#$file.6frameTranslation.fasta
#
#Can change the min orf length by changing $minORFLen
#
#######################################################################

#first get file from command prompt
my $ORIGINAL_file = shift @ARGV;
my $file;

if ($ORIGINAL_file =~ /\/common\/data\/(.+)/){
	$file = $1;
	print "$file\n";
}
else{
	$file = $ORIGINAL_file;
}


my $minORFLen = 90; #nucleotide numbering 
open OUTFASTA, ">$file.ORF.fasta";

my $outfile = "$file.6frameTranslation.fasta";
open TRANSLATION, ">$outfile";

#input file will be a fasta file of nucleotides
my $in  = Bio::SeqIO->new(-file => $ORIGINAL_file , '-format' => 'Fasta');


#loop thru each sequence in fasta file
while ( my $seq = $in->next_seq() ) {
 	my $id = $seq->id;
	my $seqLength = length($seq->seq);        
	my ($description,$translation);
   	my @frames = (0, 1, 2); 
      	
	foreach my $frame (@frames) { 
		$translation = $seq->translate(undef, undef, $frame)->seq; 
		
		#save new description
		$description = "Frame +".($frame+1);
	    
		#writting to the output file
		print TRANSLATION ">$id $description totalBP=$seqLength\n$translation\n";	
		
		
		$translation = $seq->revcom->translate(undef, undef, $frame)->seq; 
   
    		#save new description
		$description = "Frame -". ($frame+1);
	    	
		#writting to the output file
		print TRANSLATION ">$id $description totalBP=$seqLength\n$translation\n";	

  	}
  
}
close TRANSLATION;
######################################################################
#
#The second part of this program will open the file with the 6 frame translations 
#and print each ORF to a new fasta file
#
#######################################################################


#open protienFasta.txt
my $protein_in  = Bio::SeqIO->new(-file => "$outfile", '-format' => 'Fasta');

print "ProteinInFile = $outfile\n";

my $seqID;
my %longestFrame;
#loop thru each protein in this fasta

while ( my $seq = $protein_in->next_seq() ) {
	#get the id	
	$seqID = $seq->id;
  
  	#get frame and total sequence length
  	my ($frame,$seqLen);
  	my $description = $seq->description; 
  	if ($description =~ /Frame ([+|-]\d) totalBP=(\d+)/){
		$frame = $1;
		$seqLen = $2;
	}

  	#get the sequence
  	my $protein = $seq->seq; 
 
  	while ($protein =~ /\*?([\w]+\*?)/g) {
    		my $orf = $1; 
    		my $orfLen = length ($1);
    		if ($orfLen >= ($minORFLen/3) ){  
      	
			#pos returns end of m//
			my $orfStart =(pos($protein) - $orfLen);
      	
			#convert computer postion to biologically sig postion
      			$orfStart = $orfStart + 1;
    
  			my ($bp_start, $bp_stop) = aa2bp($orfStart,$orfLen, $frame, $seqLen );

			if (!defined $longestFrame{$seqID}){
				$longestFrame{$seqID} = {
					'frame' => $frame,
					'seq' => $orf,
					'length' => $orfLen,
					'bp_start' => $bp_start,
					'bp_stop' => $bp_stop,
					'seqLen' => $seqLen,
				};
			}else {
				if($longestFrame{$seqID}{length} < $orfLen){
					$longestFrame{$seqID} = {
						'frame' => $frame,	
						'seq' => $orf,
						'length' => $orfLen, 
						'bp_start' => $bp_start,
						'bp_stop' => $bp_stop,
						'seqLen' => $seqLen,
					};
				}
			}
			#print OUTFASTA ">$seqID(ORF:$bp_start..$bp_stop) Frame $frame\n$orf\n";
    		}
	}
}

#print Dumper(\%longestFrame);
foreach my $seqID(keys %longestFrame){
	my $bp_start = $longestFrame{$seqID}{bp_start};
	my $bp_stop = $longestFrame{$seqID}{bp_stop};
	my $seqLen = $longestFrame{$seqID}{seqLen};
	my $orfLen = $longestFrame{$seqID}{length} * 3;
	
	print OUTFASTA ">$seqID(ORF:$bp_start..$bp_stop) $orfLen of $seqLen  Frame $longestFrame{$seqID}{frame}  $seqID(ORF:$bp_start..$bp_stop)\n$longestFrame{$seqID}{seq}\n"; ##modified 02-16-2008
	#print OUTFASTA ">$seqID $longestFrame{$seqID}{frame}\n$longestFrame{$seqID}{seq}\n";} ## orignal print out
}
sub aa2bp {
	my ($ORFstart,$ORFlen, $frame, $seqLen) = @_;
	my $ORFstop = $ORFstart + $ORFlen - 1;
	my ($bp_start,$bp_stop);
	my $aaX3_start = $ORFstart * 3;
	my $aaX3_stop = $ORFstop * 3;
	
	if ($frame eq "+1"){
		$bp_start = $aaX3_start - 2;
		$bp_stop = $aaX3_stop + 0;
	}
	elsif ($frame eq "+2"){
		$bp_start = $aaX3_start - 1;
		$bp_stop = $aaX3_stop + 1;
	}	
	elsif ($frame eq "+3") {
		$bp_start = $aaX3_start + 0;
		$bp_stop = $aaX3_stop + 2;
	}
	elsif ($frame eq "-1"){
		$bp_stop = $seqLen - $aaX3_start + 3;
		$bp_start = $seqLen - $aaX3_stop + 1;
	}
	elsif ($frame eq "-2"){
		$bp_stop = $seqLen - $aaX3_start + 2;
		$bp_start = $seqLen - $aaX3_stop + 0;
	}
	else {
		$bp_stop = $seqLen - $aaX3_start + 1;
		$bp_start = $seqLen - $aaX3_stop - 1;
	}

	if ($bp_stop <=  0){
		$bp_stop = 1;
	}
	elsif($bp_start > $seqLen){
		$bp_start = $seqLen;
	}


	return ($bp_start, $bp_stop);
}	

