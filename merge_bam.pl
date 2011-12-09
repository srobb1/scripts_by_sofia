#!/usr/bin/perl -w
use strict;
use File::Spec;
use File::Basename;

my $dir = shift;
my $dir_path = File::Spec->rel2abs($dir);
my $lowest_dir = basename($dir_path);
print "samtools merge $dir_path/$lowest_dir.merged.bam $dir_path/*sorted.bam\n";
