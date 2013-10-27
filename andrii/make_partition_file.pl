#!/usr/bin/perl -w
use strict;

##>ccin
#(278)(451)(113)(846)(264)(266)(406)(148)(769)(304)(409)(308)(400)(231)(194)(381)(355)(268)(371)(367)(462)(92)(153)(223)(218)(227)(768)(104)
#DNA, p1=1-30
#DNA, p2=31-60

my $file =shift;
my $line = `head -n2 $file | tail -n1`;
my (@len) = $line =~ /\((\d+)\)/g;

my $total = 1;
my $i = 1;
foreach my $len (@len){
  my $start = $total;
  $total = $total + $len - 1;
  print "MOD, p$i=$start-$total\n";
  $i++; 
}
