#!/usr/bin/perl -w

use strict;

my $in_fa = shift;
my $coordinate = shift;
my $padding = shift;
if (!defined $padding or !defined $coordinate or !defined $in_fa){
  die "Please provide a fastaFile, cooridinate (ex. Chr1.28184849) and a Padding\n";
}
my ($ref , $position) = $coordinate =~ /(\S+)\.(\d+)/;
open INFA, $in_fa or die $!;
my $seq = '';
my $name;
while (my $line = <INFA>){
	chomp $line;
	if ($line =~ />(\S+)/){
		$name = $1;
		if ($seq ne '' and $name eq $ref){
			my $start =  $position-$padding-1;
			my $subseq = substr $seq, $start, ($padding*2)+1;
			print ">$name:",$start+1, "..", $position+$padding+1 ,"\n$subseq\n";
		}
	
		$seq='';	
	}else {
		$seq.=$line;
	}
}
#for last or only seq
my $start =  $position-$padding-1;
my $subseq = substr $seq, $start, ($padding*2)+1;
print ">$name:",$start+1, "..", $position+$padding+1 ," \n$subseq\n";
