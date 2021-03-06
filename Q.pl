#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#torque01.int.bioinfo.ucr.edu:
#                                                                               Req'd    Req'd      Elap
#Job ID               Username    Queue    Jobname          SessID NDS   TSK    Memory   Time   S   Time
#-------------------- ----------- -------- ---------------- ------ ----- ------ ------ -------- - --------
#4645[2639].torqu     robb        js       run_split_blast_   1628     1      4    --       --  C 00:01:58
#4645[2640].torqu     robb        js       run_split_blast_   1631     1      4    --       --  C 00:02:04
#4645[2641].torqu     robb        js       run_split_blast_   1691     1      4    --       --  C 00:02:57
#4645[2642].torqu     robb        js       run_split_blast_   2057     1      4    --       --  C 00:01:21
#4645[2643].torqu     robb        js       run_split_blast_   2118     1      4    --       --  C 00:02:40

my $now_string = localtime;
print "$now_string:\n";
my @jobs = `qstat -urobb -t`;
my %status;
foreach  my $line (@jobs){
  chomp $line;
  next if $line !~ /^\d/;
  my @line = split /\s+/ , $line;
  my $proc = $line[6] ne '--' ? $line[6] : 1;
  $status{$line[2]}{$line[9]}+=$proc;
}
print Dumper \%status;
