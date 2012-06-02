#!/usr/bin/perl
# ejr - 20120216
# open two fastq files, write matching pairs to right and lefts files and then
# unpaired to a third file
# inefficient IO, but has a lower memory footprint for big files
# I should find a better way to do this. BerkeleyDb and SQlite are slower
use strict;
use warnings;


my $rightfile = $ARGV[0] . "_rights";
my $leftfile = $ARGV[1] . "_lefts";
my $unpairedfile = $ARGV[0] . "_unpaired";
my $dbfile = $ARGV[0] . "_db";
my %lefts;
my %rights;

open (LEFT, $ARGV[0]) or die "cannot open file:$!\n";
open (RIGHT, $ARGV[1]) or die "cannot open file:$!\n";
open (LEFTOUT, ">$rightfile") or die "cannot open file:$!\n";
open (RIGHTOUT, ">$leftfile") or die "cannot open file:$!\n";
open (UNPAIRED, ">$unpairedfile") or die "cannot open file:$!\n";

while (my $left_header = <LEFT>) {
       $lefts{$left_header} = 1;
       <LEFT>;
       <LEFT>;
       <LEFT>;
}

close LEFT;

while (my $right_header = <RIGHT>) {
       my $right_seq = <RIGHT>;
       my $right_qual_header = <RIGHT>;
       my $right_qual = <RIGHT>;
       my $right_full = $right_header.$right_seq.$right_qual_header.$right_qual;
           if (exists($lefts{$right_header})){
               print RIGHTOUT $right_full;
               $rights{$right_header} = 1;
           } else {
               print UNPAIRED $right_full;
           }
}
close RIGHT;

open (LEFT, $ARGV[0]) or die "cannot open file:$!\n";

while (my $left_header = <LEFT>) {
       my $left_seq = <LEFT>;
       my $left_qual_header = <LEFT>;
       my $left_qual = <LEFT>;
       my $left_full =$left_header.$left_seq.$left_qual_header.$left_qual;
       if ($rights{$left_header}) {
           print LEFTOUT $left_full;
       } else {
           print UNPAIRED $left_full;
       }
}
