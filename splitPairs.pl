#!/usr/bin/perl -w
# Liping Hou

use strict;

die("Usage split_fastq.pl out.prefix in.fastq\n") if (@ARGV == 0);
my $pre = shift(@ARGV);
my ($fq1, $fq2);
open($fq1, ">${pre}_1.fastq") || die;
open($fq2, ">${pre}_2.fastq") || die;
while (<>) {
    chomp;
    print $fq1 $_ . "\n";
    $_ = <>; print $fq1 $_;
    $_ = <>; print $fq1 $_;
    $_ = <>; print $fq1 $_;
    $_ = <>; chomp; print $fq2 $_ . "\n";
    $_ = <>; print $fq2 $_;
    $_ = <>; print $fq2 $_;
    $_ = <>; print $fq2 $_;
}
close($fq1);
close($fq2);

