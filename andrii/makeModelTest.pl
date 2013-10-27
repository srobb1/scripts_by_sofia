#!/usr/bin/perl -w
use strict;

my $dir = shift;

my @files = <$dir/*phy>;
foreach my $file (@files){
  my @file = split /\// , $file;
  my $name = pop @file;
  open SH, ">$file.modeltest.sh";


print SH"
\#This is the modeltest 
\#PBS -l nodes=1:ppn=8 -q highmem

phy=$name

module load RAxML/7.3.2
cd /rhome/andriig/Phylo
mkdir -p ModelTest
cd ModelTest
ln -s ../\$phy
ProteinModelSelection.pl \$phy >& $name.modelselection.out
";

}
