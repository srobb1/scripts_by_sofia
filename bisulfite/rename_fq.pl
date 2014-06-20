#!/usr/bin/perl -w
use strict;
use File::Copy;

my $dir = shift;  #dir of fastq files
my $file = shift; #rename file# old_name  new name

if (!defined $dir or !defined $file){
  die "Please provide the 
	1. directory of fastq files to be renamed and 
	2. a file with \"old_name new_name\

    and a symlink of your file will be generated in your current working directory\n";
 
}
open IN, $file or die "Can't open Rename file for reading: $file\n";
my %rename;
while (my $line = <IN>){
  chomp $line;
  my ($old,$new) = split /\s+/, $line;
  $rename{$old}=$new;
}


my @fastqs = <$dir/*.fastq>;
foreach my $fastq (@fastqs){
  my @path = split /\// , $fastq;
  my $old = pop @path; 
  my $new = $rename{$old};
  symlink $old,$new;
}
