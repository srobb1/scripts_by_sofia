#!/usr/bin/perl -w

use strict;
my $file = shift;
open (INFILE, $file) or die "Can't open $file $!\n";

=cut
while (my $line = <INFILE>){
        chomp $line;
        my @line = split /\t/, $line;
        my $left_seq = $line[4];
        my $right_seq = $line[5];
        $left_seq =~ s/left_flanking_seq=//;
        $right_seq =~ s/right_flanking_seq=//;
        print ">A119.$line[0].$line[1]\n$left_seq$right_seq\n";
}
=cut

while (my $line = <INFILE>){
	chomp $line;
	my @line = split /\t/, $line;
	my $left_seq = $line[5];
	my $right_seq = $line[6];
	$left_seq =~ s/left_flanking_seq=//;
	$right_seq =~ s/right_flanking_seq=//;
	print ">EG4.$line[1].$line[2]\n$left_seq$right_seq\n";
}
