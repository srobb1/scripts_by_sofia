#!/usr/bin/perl -w
use strict;
use File::Spec;

my $dir         = shift;
die "Please provide a directory of bam files" if !defined $dir and !-d $dir;
$dir =~ s/\/$//;

my $bam_dir    = File::Spec->rel2abs($dir);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);


my @bam_files = <$bam_dir/*bam>;

foreach my $bam (@bam_files){
  my ($base) = $bam =~ /(\S+)\.bam/;
  my @base = split /\// , $base;
  my $name = pop @base;
  open SH , ">$name.getUniqReads.sh"; 
  print SH "samtools view -q 5 $bam > $base.uniq.sam\n";
}


##now need to run sam2fq on the directory of sam files.
