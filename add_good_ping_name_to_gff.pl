#!/usr/bin/perl -w
use strict;
use Data::Dumper;
my $file = shift;
my @pings= qw (
ID=HEG4.ping.te_insertion_site.Chr1.4220012;Name=Ping2
ID=EG4.ping.te_insertion_site.Chr1.4220012;Name=Ping2
ID=HEG4.ping.te_insertion_site.Chr3.28019802;Name=Ping3
ID=EG4.ping.te_insertion_site.Chr3.28019802;Name=Ping3
ID=HEG4.ping.te_insertion_site.Chr7.26460309;Name=Ping4
ID=EG4.ping.te_insertion_site.Chr7.26460309;Name=Ping4
ID=HEG4.ping.te_insertion_site.Chr9.10863120;Name=Ping5
ID=EG4.ping.te_insertion_site.Chr9.10863120;Name=Ping5
ID=HEG4.ping.te_insertion_site.Chr9.16690614;Name=Ping7
ID=EG4.ping.te_insertion_site.Chr9.16690614;Name=Ping7
ID=HEG4.ping.te_insertion_site.Chr9.13736143;Name=Ping6
ID=EG4.ping.te_insertion_site.Chr9.13736143;Name=Ping6
ID=HEG4.ping.te_insertion_site.Chr1.2640502;Name=Ping1
ID=EG4.ping.te_insertion_site.Chr1.2640502;Name=Ping1
ID=A119.ping.te_insertion_site.Chr1.2640502;Name=Ping1
ID=A123.ping.te_insertion_site.Chr1.2640502;Name=Ping1
ID=A119.ping.te_insertion_site.Chr1.31034126;Name=Ping9
ID=A119.ping.te_insertion_site.Chr2.27861999;Name=Ping10
ID=A119.ping.te_insertion_site.Chr2.5763675;Name=Ping11
ID=A119.ping.te_insertion_site.Chr3.15950469;Name=Ping12
ID=A119.ping.te_insertion_site.Chr4.5816107;Name=Ping13
ID=A119.ping.te_insertion_site.Chr5.28565219;Name=Ping14
ID=A123.ping.te_insertion_site.Chr1.33282259;Name=Ping15
ID=A123.ping.te_insertion_site.Chr10.1267746;Name=Ping16
ID=A123.ping.te_insertion_site.Chr11.1807322;Name=Ping17
ID=A123.ping.te_insertion_site.Chr2.12298998;Name=Ping18
ID=A123.ping.te_insertion_site.Chr2.30902247;Name=Ping19
ID=A123.ping.te_insertion_site.Chr2.31053564;Name=Ping20
ID=A123.ping.te_insertion_site.Chr5.171255;Name=Ping21
ID=A123.ping.te_insertion_site.Chr5.334874;Name=Ping22
ID=A123.ping.te_insertion_site.Chr8.24881082;Name=Ping23
ID=NB.ping.te_insertion_site.Chr6.23526981;Name=Ping8 );
my %rename;
foreach my $ping (@pings){
  #print $ping,"\n";
  my ($loc , $name) = $ping =~ /(Chr\d+\.\d+);Name=(Ping\d+)/;
  $rename{$loc}=$name;
}
#print Dumper \%rename;
open IN, $file or die "Can't open $file\n";
while (my $line = <IN>){
  chomp $line;
  #Name=ping.Chr7_26460307_26460309;Note
  if ($line =~ /Name=(ping\.(Chr\d+)_\d+_(\d+));/){
   #print "$1:: $2.$3\n";
   my $new = $rename{"$2.$3"};
   #print "new:$new\n";
   $line =~ s/Note=/Note=$new,/;
   print $line,"\n";
  }
}
