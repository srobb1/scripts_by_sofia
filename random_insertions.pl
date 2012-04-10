#!/usr/bin/perl -w
use strict;

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

foreach my $chr (sort keys %genome){
  my $size = $genome{$chr};
  for (my $i = 0 ; $i < 2000 ; $i++) {
    my $insert = get_random_bp($size);
     print 
"$chr\tcontrol\ttransposable_element_attribute\t$insert\t$insert\t.\t.\t.\tID=$chr.$insert;avg_flankers=0;spanners=0;type=0;TE=control;TSD=control\n";
  }
}


sub get_random_bp {
  my $max = shift;
  return int(rand($max))
}
