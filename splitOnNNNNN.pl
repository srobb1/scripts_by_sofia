#!/usr/bin/perl
use strict;
$/='>';

my $fa_file = shift;
my $pattern = 'N{1000,}';
my $len = 3;
my $insert_bp = 3; # following bp-6 (A) before bp-7 the mping will insert

open INFA , $fa_file or die "Can't open $fa_file\n";
<INFA>; #throw out first '>'
while (my $line = <INFA>){
   chomp $line;
   my ($header,@seq) = split /\n/ , $line;
   my ($id , $desc) = $header =~ /(\S+)\s*(.*)/;
   my $seq = join '' , @seq;
   my $count; 
   while ($seq =~ /$pattern/gi){
     $count++
   }
   my @subseqs = split /$pattern/ , $seq;
   my $i = 1;
   foreach my $subseq (@subseqs){
     if (length $subseq > 0 ){
       print ">$id-$i\n$subseq\n";
       $i++;
     }
  }
}
