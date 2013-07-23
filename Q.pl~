#!/usr/bin/perl -w
use strict;
use Data::Dumper;
## echo "running:" ; qstat -u robb -t | awk {'print $10'} | grep -c R ;  echo "Waiting:" ; qstat -urobb -t | awk {'print $10'}| grep -c Q ; echo "Holding:" ; qstat -urobb -t | awk {'print $10'} | grep -c H 


my $now_string = localtime;
print "$now_string:\n";
my @jobs = `qstat -urobb -t`;
my %status;
foreach  my $line (@jobs){
  chomp $line;
  next if $line !~ /^\d/;
  my @line = split /\s+/ , $line;
  $status{$line[2]}{$line[9]}++;
}
print Dumper \%status;
