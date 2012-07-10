#!/usr/bin/perl -w

use strict;

die("Usage split_fastq.pl out.prefix in.fastq\n") if (@ARGV == 0);
my $pre = shift(@ARGV);
my ($fq1, $fq2);
open($fq1, ">${pre}_1.fastq") || die;
open($fq2, ">${pre}_2.fastq") || die;
while (my $header = <>) {
   unless ($header =~ /^@/) {
      die(
"doesn't look like a header, got line:\n$header  but expected line to start with '\@'\n"
      );
    }
    print $fq1 $header;
    for ( 0 .. 2 ) {
      my $line = <>;
      print $fq1 $line;
    }
   $header         = <>; 
   unless ($header =~ /^@/) {
      die(
"doesn't look like a header, got line:\n$header  but expected line to start with '\@'\n"
      );
    }
    print $fq2 $header;
    for ( 0 .. 2 ) {
      my $line = <>;
      print $fq2 $line;
    }
}
close($fq1);
close($fq2);
