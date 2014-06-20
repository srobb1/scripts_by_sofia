#!/usr/bin/perl -w
use strict;

#mping.Chr10_22433383_22433385 mping.fwd
#mping.Chr10_20589625_20589627 mping.fwd
#mping.Chr10_17020823_17020825 mping.fwd
#mping.Chr8_23643362_23643364 mping.rc

my $file = 'mpingRef.list';
open IN, $file or die "Can't open $file\n";
while (my $line = <IN>){
  chomp $line;
  my ($ref,$dir) = split /\s/ , $line;
  my $strand = $dir =~ /fwd/ ? '+' : '-';
  
  print join("\t",$ref,"spring2014","region",1,1630,'.','.','.',"ID=$ref;name=$ref"),"\n";
  print join("\t",$ref,"spring2014","transposable_element",601,1030,'.',$strand,'.',"ID=$ref.601..1030;name=mping"),"\n";
}
