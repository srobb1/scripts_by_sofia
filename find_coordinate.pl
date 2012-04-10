#!/usr/bin/perl -w
use strict;

$/ = '>';
my $file = shift;
my $pattern = shift;

open INFASTA , $file;
<INFASTA>;
while (my $line = <INFASTA>){
  chomp $line;
  my ($header , @seq) = split /\n/ , $line;
  my $seq = join '' , @seq;
  my ($id) = $header =~ /^(\S+)/;

  if ($seq =~ /$pattern/g){
    my $end = pos ($seq);
    my $pattern_len = length $pattern;
    my $start = $end - $pattern_len + 1;
    print "$id($pattern):$start..$end\n";
  }else {
    #print "$id($pattern):NOTFOUND\n";
  }


}

