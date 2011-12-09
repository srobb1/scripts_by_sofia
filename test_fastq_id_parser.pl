#!/usr/bin/perl -w
use strict;

my $fq_file= shift;

open (FASTQ, '<', $fq_file) or die $!;

my $line_count = 4;
while (<FASTQ>){
	#print "$line_count: $_\n";
	chomp;
        if (/^@/ and $line_count == 4) {
		print "$_\n";
        	#print "$_\n";
		$line_count = 0;
	}
    	$line_count++;	
}
