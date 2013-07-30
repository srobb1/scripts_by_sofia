#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#torque01.int.bioinfo.ucr.edu:
#                                                                               Req'd    Req'd      Elap
#Job ID               Username    Queue    Jobname          SessID NDS   TSK    Memory   Time   S   Time
#-------------------- ----------- -------- ---------------- ------ ----- ------ ------ -------- - --------
#4295.torque01.in     xzhang      batch    Arfima_11_04_1_n    --    --     --     --  168:00:0 R 163:16:2   n28/3
#4296.torque01.in     xzhang      batch    Arfima_11_04_2_n    --    --     --     --  168:00:0 R 163:10:4   n29/1
#4297.torque01.in     xzhang      batch    Arfima_32_04_1_n    --    --     --     --  168:00:0 R 163:11:2   n29/3
#4299.torque01.in     xzhang      batch    Arfima_11_04_1_p    --    --     --     --  168:00:0 R 163:10:2   n29/5
#4300.torque01.in     xzhang      batch    Arfima_11_04_2_p    --    --     --     --  168:00:0 R 163:11:0   n28/2
#4336.torque01.in     xlu006      batch    tss.sh              --    --     --     --  168:00:0 R 146:38:4   n33/30
#4486.torque01.in     croberts    batch    STDIN               --      1     40    --  168:00:0 R 119:03:1   n33/47+n33/46+n33/45+n33/44+n33/43+n33/42+n33/41+n33/40+n33/39+n33/38+n33/37+n33/36+n33/35+n33/34+n33/33+n33/32+n33/31+n33/28+n33/27+n33/26+n33/25+n33/24+n33/23+n33/22+n33/21+n33/20+n33/19+n33/18+n33/17+n33/16+n33/15+n33/14+n33/13+n33/12+n33/11+n33/10+n33/9+n33/8+n33/2+n33/1

my $now_string = localtime;
print "$now_string:\n";
my @jobs = `qstat -n1`;
my %status;
foreach  my $line (@jobs){
  chomp $line;
  next if $line !~ /^\d/;
  my @line = split /\s+/ , $line;
  my @nodes = split /\+/ , $line[11];
  foreach my $node (@nodes){
    $node =~ s/\/\d+//;
    $status{$line[2]}{$line[9]}{$node}++;
  }
  #my $proc = $line[6] ne '--' ? $line[6] : 1;
  #$status{$line[2]}{$line[9]}+=$proc;
}
print Dumper \%status;

__END__
foreach my $q (sort keys %status){
  print "$q =>\n";
  foreach my $state (sort keys %{$status{$q}}){
    print "\t$state =>\n";
    foreach my $node (sort keys %{$status{$q}{$state}}){
      my $count = $status{$q}{$state}{$node};
      if ($state ne 'Q'){
        print "\t\t$node\t$count\n";
      }else{
        print "\t\t\t$count\n";
      }
    }
  }
}
