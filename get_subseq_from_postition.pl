#!/usr/bin/perl -w

use strict;

my $in_fa = shift;
my $position = shift;
my $padding = 300;
open INFA, $in_fa or die $!;
my $seq = '';
my $name;
while (my $line = <INFA>){
	chomp $line;
	if ($line =~ />(\S+)/){
		$name = $1;
		if ($seq ne ''){
			my $start =  $position-$padding;
			my $subseq = substr $seq, $start, 600;
			print ">$name\n$subseq\n";
		}
	
		$seq='';	
	}else {
		$seq.=$line;
	}
}
#for last or only seq
my $start =  $position-$padding;
my $subseq = substr $seq, $start, 600;
print ">$name\n$subseq\n";
