#!/usr/bin/perl -w
use strict;

my $insert_count = shift;


##MSUr7
my %genome = (
  'Chr1' => 43270923 ,
  'Chr2' => 35937250 ,
  'Chr3' => 36413819 ,
  'Chr4' => 35502694 ,
  'Chr5' => 29958434, 
  'Chr6' => 31248787 ,
  'Chr7' => 29697621 ,
  'Chr8' => 28443022 ,
  'Chr9' => 23012720 ,
  'Chr10' => 23207287 ,
  'Chr11' => 29021106 ,
  'Chr12' => 27531856
);

##generate 1000 insertions per chromosome


for (my $i = 0 ; $i < $insert_count ; $i++){
  my $chr_num = get_random_int(12) + 1;
  my $size = $genome{"Chr$chr_num"};
  my $insert = get_random_int($size) + 1;
  print "Chr$chr_num:$insert\n"; 
}


sub get_random_int {
  my $max = shift;
  return int(rand($max));
}
