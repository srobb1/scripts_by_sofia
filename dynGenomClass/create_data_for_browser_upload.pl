#!/usr/bin/perl -w

use strict;
my $file = shift;
##
open INFILE, $file or die $!;
<INFILE>;

my $last_chr = ''; 
while (my $line = <INFILE>){
        chomp $line;
        my ($source, $chr, $pos) = split /\t/ ,$line;
        if ($chr ne $last_chr){
           open OUTFILE, ">$source.$chr.forBrowser.txt";
	   print OUTFILE "
[A123_insertion]
glyph       = triangle
citation    = mPing Insertion in A123
key         = mPing Insertion in A123

";
        }
        print OUTFILE $source."_insertion\t$chr:$pos\t$chr:$pos..$pos\n";
        $last_chr = $chr;
}
