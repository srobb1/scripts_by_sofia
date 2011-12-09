#!/usr/bin/perl -w

use strict;
use File::Spec;

#cd cromosome01
#mkdir unpaired ; mv *unPaired* unpaired
#cat *1.fq > A157_chromosome01_combined_1.fq ; cat *2.fq > A157_chromosome01_combined_2.fq ; mv A157* ../.
#cd unpaired ; cat *fq > A157_chromosome01_combined_unPaired_1.fq ; mv *combined* ../../. ; cd ../../
my $lib = shift;
my $chromosome = shift;
my $dir = shift;

my $dir_path = File::Spec->rel2abs($dir);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);
open OUTSH, ">$dir_path/$chromosome.combine_fq.sh";
print OUTSH "#!/bin/bash\n\n";

print OUTSH "cd $dir_path\n";
print OUTSH "mkdir $dir_path/unpaired\n";
print OUTSH "mv $dir_path/p*unPaired*.fq $dir_path/unpaired/.\n";
print OUTSH "cat $dir_path/p*1.fq > $dir_path/$lib"."_".$chromosome."_combined_1.fq\n";
print OUTSH "cat $dir_path/p*2.fq > $dir_path/$lib"."_".$chromosome."_combined_2.fq\n";
print OUTSH "cd $dir_path/unpaired\n";
print OUTSH "cat $dir_path/unpaired/*fq > $dir_path/$lib"."_".$chromosome."_combined_unPaired.fq\n";
