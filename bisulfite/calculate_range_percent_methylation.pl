#!/usr/bin/perl -w
use strict;
my $type   = shift;    ## CHH, CHG, CpG
my $strain = shift;
my $range  = shift;
my $dir    = shift;##split_methylation

my ( $ref, $range_start, $range_end ) = $range =~ /(.+)\:(\d+)\.\.(\d+)/;
my %methy;
open IN, "$dir/${strain}_${ref}_${type}.split.txt";    #NB_Chr9_CpG.split.txt
while ( my $line = <IN> ) {
  chomp $line;
  my ( $pos, $methyl_yes, $methyl_no ) = split /\t/, $line;
  next unless $pos >= $range_start;
  next unless $pos <= $range_end;
  last if $pos > $range_end;
  $methyl_yes = defined $methyl_yes ? $methyl_yes : 0;
  $methyl_no  = defined $methyl_no  ? $methyl_no  : 0;
  $methy{$pos}{$type}{yes} += $methyl_yes;
  $methy{$pos}{$type}{no}  += $methyl_no;
}
my $region_yes = 0;
my $region_no  = 0;
my $to_print = "strain\ttype\tref\tpos\t%\ttotal\tmethy_yes\tmethy_no\n";
foreach my $pos ( sort { $a <=> $b } keys %methy ) {
  foreach my $type ( keys %{ $methy{$pos} } ) {
    my ( $methy_yes, $methy_no ) = ( 0, 0 );
    if ( exists $methy{$pos}{$type}{yes} ) {
      $methy_yes = $methy{$pos}{$type}{yes};
      $region_yes += $methy_yes;
    }
    if ( exists $methy{$pos}{$type}{no} ) {
      $methy_no = $methy{$pos}{$type}{no};
      $region_no += $methy_no;
    }
    my $total = $methy_yes + $methy_no;
    my $percent = ( $methy_yes / ($total) ) * 100;

    #if ( $percent == 0 ) {
    #  $percent = -100;
    #}
    my $pretty_percent = sprintf( '%.1f', $percent );
    $to_print .=
      "$strain\t$type\t$ref\t$pos\t$pretty_percent\t$total\t$methy_yes\t$methy_no\n";
  }
}
my $region_total          = $region_yes + $region_no;
my $region_percent        = ( $region_yes / $region_total ) * 100;
my $pretty_region_percent = sprintf( '%.1f', $region_percent );

print "$range $pretty_region_percent% $type in $strain\n";
print "=========================================================\n";
print $to_print;
print "=========================================================\n";
print "$range $pretty_region_percent% $type in $strain\n";
