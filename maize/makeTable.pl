#!/usr/bin/perl -w
use strict;
##to use sites
my $file = shift;
my @NAMs;
while (my $NAM = <DATA>){
  chomp $NAM;
  next if $NAM =~ /^\s+$/;
  push @NAMs , $NAM;
}
@NAMs = sort @NAMs;

open IN , $file or die "Can't open $file\n";
my %sites;
while (my $line = <IN>){
  chomp $line;
  next if $line !~ /\d+\|/; 
  my ($NAM, $site) = split /\t/ , $line;
  $sites{$site}{$NAM}=1; 
}

print "site\t" , join ("\t" , @NAMs) , "\n";
foreach my $site (sort keys %sites ){
  my @found;
  foreach my $NAM (@NAMs){
    if (exists $sites{$site}{$NAM}){
      push @found , 1;
    }else {
      push @found , 0;
    }
  } 
    print "$site\t" , join ("\t" , @found) , "\n";
}

__DATA__
CML247
Ki3
Tzi8
Il14H
Ky21
CML52
P39
B97
NC358
M162W
B73
