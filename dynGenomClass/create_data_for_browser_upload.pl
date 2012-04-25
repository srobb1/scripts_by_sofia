#!/usr/bin/perl -w

## use the class.query.out file: 
## source  ref     start   insert_type     spanners        avg.flankers    insert_feature  f_start f_end   f_strand        f_name  f_note
## RIL39_2 Chr10   11955072        homozygous      0       21.5    intergenic      11951354        11961940        LOC_Os10g22960:11947165..11951353(+)/LOC_Os10g22970:11961941..11962177(-)       n/a     n/a
use strict;
my $file = shift;
open INFILE, $file or die $!;
<INFILE>;

my $last_chr = ''; 
while (my $line = <INFILE>){
        chomp $line;
        my ($source, $chr, $pos, $class, $spanners, $avg_flank) = split /\t/ ,$line;
        if ($chr ne $last_chr){
           open OUTFILE, ">$source.$chr.forBrowser.txt";
	   print OUTFILE "
[$source.insertion]
glyph       = triangle
citation    = mPing Insertion in $source
key         = mPing Insertion in $source

";
        }
        print OUTFILE $source."_insertion\t$chr:$pos\t$chr:$pos..$pos\tNote=$class,spanners($spanners),avg.flankers($avg_flank)\n";
        $last_chr = $chr;
}
