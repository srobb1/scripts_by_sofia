#!/usr/bin/perl -w
use strict;

my @files = @ARGV;

my %comparison;

##TE      TSD     Exper   chromosome      insertion_site  left_flanking_read_count        right_flanking_read_count       left_flanking_seq       right_
foreach my $file (@files){
  open IN , "$file" or die "Can't open $file\n";
  while (my $line = <IN>){
    chomp $line;
    next if $line =~ /Exper/;
    my ($TE, $tsd, $strain, $chr, $pos, $left_count, $right_count, @rest) = split /\t/, $line;
    push @{$comparison{$TE}{"$chr.$pos"}} , $strain;
  }
}
my %counts;
print "TE\tloc\tcount_of_libraries_sharing_this_insert\n";
foreach my $TE ( sort keys %comparison){
  foreach my $loc (sort keys %{$comparison{$TE}}){
    my $strains_str = join (',', sort @{$comparison{$TE}{$loc}});
    $counts{$TE}{$strains_str}++;
    my $count_shared = scalar @{$comparison{$TE}{$loc}};
    print "$TE\t$loc\t$count_shared\t$strains_str\n";
  }
}
print "\n\ncounts\n";
print "TE\tstrain\tuniqCount\n";
foreach my $TE (sort keys %counts){
  foreach my $strains (sort keys %{$counts{$TE}}){
    print "$TE\t$strains\t$counts{$TE}{$strains}\n";
  }
}
