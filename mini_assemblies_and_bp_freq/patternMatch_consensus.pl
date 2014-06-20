#!/usr/bin/perl
use strict;

## find all sites in genome that have the consensus mping insertion sequence.
## pattern is 3bp in length
## insertion bp is 3


$/='>';

my $fa_file = shift;
my $pattern = 'T[TA]A';
my $len = 3;
my $insert_bp = 3; # following bp-6 (A) before bp-7 the mping will insert

open INFA , $fa_file or die "Can't open $fa_file\n";
<INFA>; #throw out first '>'
while (my $line = <INFA>){
   chomp $line;
   my ($header,@seq) = split /\n/ , $line;
   my ($id , $desc) = $header =~ /(\S+)\s*(.*)/;
   my $seq = join '' , @seq;
   
   while ($seq =~ /($pattern)/g){
     my $pos = pos $seq;
     my $insert_loc = $pos - ($len - $insert_bp);
     print "$id:$insert_loc\n";

   }
}
