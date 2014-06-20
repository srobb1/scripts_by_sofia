#!/usr/bin/perl -w

use strict;


my $file = shift;

open INFILE, $file or die $!; 
open OUTFILE, ">$file.uniq" or die $!;

my %seqIDs;

while (my $seqID = <INFILE>){
	my $seq = <INFILE>;
	my $qualHeader = <INFILE>;
	my $qual = <INFILE>;
	if (!exists $seqIDs{$seqID}){
		$seqIDs{$seqID} = 1;
		print OUTFILE $seqID ,  $seq , $qualHeader, $qual;
	}
	#else {
		#don't print to outfile and move on to next seq record
	#}
}
