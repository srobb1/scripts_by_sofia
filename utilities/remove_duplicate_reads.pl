#!/usr/bin/perl -w
use strict;

# takes a fq file
# writes to STDOUT
# writes out a file in the same order minus the duplicates
# so.. the second or more copy will be absent in the output

my $file = shift;
my %dupCheck;
open INFILE, $file or die "Can't open $file $!\n";

while (my $header = <INFILE>){
 	chomp $header;
	my $seq = <INFILE>;
	chomp $seq;
	my $qual_header = <INFILE>;
	chomp $qual_header;
	my $qual = <INFILE>;
	chomp $qual;
	next if exists $dupCheck{$header}{$seq};
        $dupCheck{$header}{$seq} = 1;
	print "$header\n$seq\n$qual_header\n$qual\n";
}
