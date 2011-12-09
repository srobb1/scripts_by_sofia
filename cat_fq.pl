#!/usr/bin/perl -w
use strict;
use File::Spec;
#combines all the fq files in one directory
#provide the direcotry name
## for multiple directories try this:
## for i in `ls` ; do cat_fq.pl $i ; done
my $dir = shift;
my $dir_path = File::Spec->rel2abs($dir);

my @dirs = split '/' , $dir_path;
my $lowest_dir = pop @dirs;
my $one_up = join '/' , @dirs;

my $cat1 = "cat $dir_path/*_1.fq > $one_up/$lowest_dir"."_1.fq" ; 
`$cat1` ; 
my $cat2 = "cat $dir_path/*_2.fq > $one_up/$lowest_dir"."_2.fq"; 
`$cat2`; 
my $cat3 = "cat $dir_path/*unPaired.fq > $one_up/$lowest_dir"."_unPaired.fq";
`$cat3`;
