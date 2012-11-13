#!/usr/bin/perl 
use strict;
use warnings;
use Data::Dumper;

my %depth;
my $file = shift;
my $db   = shift;
my %peaks;
print "group_size_range\tgroup_size\tmean\tgroup_mean\tgroup_max\tpos\n";
open INFILE, $file;
my $last_chr = 'Chr1';

while ( my $line = <INFILE> ) {
  chomp $line;
  my ( $chr, $pos, $depth ) = split /\s+/, $line;
  push @{ $depth{$chr}{depth} }, $depth;
  if ( $chr ne $last_chr ) {
    my $sum;
    $sum += $_ for @{ $depth{$last_chr}{depth} };
    my $mean = $sum / @{ $depth{$last_chr}{depth} };
    $depth{$last_chr}{mean} = $mean;

    findPeaks($last_chr);
    delete $depth{$last_chr};
    $last_chr = $chr;
  }
}

sub findPeaks {
  my $chr          = shift;
  my $group        = 0;
  my $gt_mean_prev = '';
  my @mins_maxs;
  my $i = 0;
  my %peaks;
  my %groups;
  my $mean = $depth{$chr}{mean};
  for my $d ( @{ $depth{$chr}{depth} } ) {
    $i++;
    my $gt_mean = $d > $mean ? 1 : 0;

    unless ( $gt_mean eq $gt_mean_prev ) {
      $gt_mean_prev = $gt_mean;
      $group++;
      $mins_maxs[$group] = $d;
    }

    if ($gt_mean) {
      $mins_maxs[$group] = $d if $d > $mins_maxs[$group];
    }
    else {
      $mins_maxs[$group] = $d if $d < $mins_maxs[$group];
    }
    push @{ $peaks{$group}{pos} },    $i;
    push @{ $peaks{$group}{values} }, $d;
    $peaks{$group}{chr}           = $chr;
    $peaks{$group}{group_mean}    = $gt_mean;
    $peaks{$group}{group_min_max} = $mins_maxs[$group];
  }

  foreach my $group ( sort { $a <=> $b } keys %peaks ) {
    my $group_size  = @{ $peaks{$group}{pos} };
    my $chr         = $peaks{$group}{chr};
    my $mean        = $depth{$chr}{mean};
    my $num_len     = length $group_size;
    my $first_digit = substr $group_size, 0, 1;
    my $bin         = $first_digit . '0' x ( $num_len - 1 );
    if ( $group_size < 100 ) {
      $bin = 0;
    }
    my $group_sum = 0;
    foreach my $val ( @{ $peaks{$group}{values} } ) {
      $group_sum += $val;
    }
    my $group_mean    = $group_sum / $group_size;
    my $group_min_max = $peaks{$group}{group_min_max};
    my $group_start   = ${ $peaks{$group}{pos} }[0];
    my $group_end     = ${ $peaks{$group}{pos} }[-1];
    if ( $group_size > 10 and $group_mean > $mean + 200 ) {
      $groups{$bin}{"$chr:$group_start..$group_end"}{chr}  = $chr;
      $groups{$bin}{"$chr:$group_start..$group_end"}{mean} = $group_mean;
      $groups{$bin}{"$chr:$group_start..$group_end"}{max}  = $group_min_max;
      $groups{$bin}{"$chr:$group_start..$group_end"}{size} = $group_size;
    }
  }
  my @ranges;
  foreach my $bin ( sort { $a <=> $b } keys %groups ) {
    foreach my $pos (
      sort { $groups{$bin}{$a}{mean} <=> $groups{$bin}{$b}{mean} }
      keys %{ $groups{$bin} }
      )
    {

      push @ranges, $pos;
      my $chr    = $groups{$bin}{$pos}{chr};
      my $mean   = $depth{$chr}{mean};
      my $end    = $bin + 99;
      my $g_mean = $groups{$bin}{$pos}{mean};
      my $size   = $groups{$bin}{$pos}{size};
      my $g_max  = $groups{$bin}{$pos}{max};
      print "$bin-$end\t$size\t$mean\t$g_mean\t$g_max\t$pos\n";
    }

  }

  my $ranges = join( ' ', @ranges );
  warn
"/home_stajichlab/robb/bin/get_subseq_from_range.pl $db $ranges > peaks.fa\n";
}
