#!/usr/bin/perl -w
use strict;


my $file = shift;
##>primersetName
##PRIMERSEQ1nnnnnnnnnnnnPRIMERSEQ2

`blastall -p blastn -i $file -d  ~/Wessler-Rice/Genome/index/MSU_r7.all.fa -o primer.blastout -m8`;


my %primers;
my $last = '';
my $count;
open BLASTOUT, "primer.blastout" or die "Can't open primer.blastout";
#1	Chr11	100.00	21	0	0	45	65	21964649	21964629	0.004	42.1
#1	Chr11	100.00	16	0	0	1	16	21963880	21963895	3.6	32.2
while ( my $line = <BLASTOUT> ) {
  chomp $line;
  my (
    $primer_pair, $subject,     $perId,  $alignLen,
    $mm,          $gapOpenings, $qStart, $qEnd,
    $sStart,      $sEnd,        $e,      $score
  ) = split /\t/, $line;
  if ( $primer_pair ne $last ) {
    $count = 0;
  }
  else {
    $count++;
  }
  $primers{$primer_pair}{$count}{sub}    = $subject;
  $primers{$primer_pair}{$count}{alnLen} = $alignLen;
  $primers{$primer_pair}{$count}{mm}     = $mm;
  $primers{$primer_pair}{$count}{perID}  = $perId;
  $primers{$primer_pair}{$count}{qStart} = $qStart;
  $primers{$primer_pair}{$count}{qEnd}   = $qEnd;
  $primers{$primer_pair}{$count}{sStart} = $sStart;
  $primers{$primer_pair}{$count}{sEnd}   = $sEnd;
  $primers{$primer_pair}{$count}{e}      = $e;
  $primers{$primer_pair}{$count}{line}   = $line;
  $last                                  = $primer_pair;
}
my %warn;
my %product_size;
foreach my $primer_set ( keys %primers ) {
  my $count = scalar keys %{ $primers{$primer_set} };
  if ( $count < 2 ) {
    push @{ $warn{$primer_set} } , 'only 1 primer mapped';
  }
  elsif ( $count > 2 ) {
    push @{ $warn{$primer_set} } ,
      'more than 1 hit for one or more of your primers in this set';
  }
  elsif ( $count == 2 ) {
    my $sub_1 = '';
    my @values;
    foreach my $count ( keys %{ $primers{$primer_set} } ) {
      my $sub = $primers{$primer_set}{$count}{sub};
      if ( $sub_1 eq '' ) {
        $sub_1 = $sub;
        push @values, $primers{$primer_set}{$count}{sStart} , $primers{$primer_set}{$count}{sEnd};
      }
      elsif ( $sub eq $sub_1 ) {
        push @values, $primers{$primer_set}{$count}{sStart} , $primers{$primer_set}{$count}{sEnd};
      }
    }
    my @sorted_values = sort {$a <=> $b} @values;
    my $smallest = shift @sorted_values;
    my $largest = pop @sorted_values;
    if (!exists $warn{$primer_set}){
      print "$primer_set\t",$largest-$smallest+1,"\n";
    }else{
      print "$primer_set\t",join (",",@{$warn{$primer_set}}),"\n"; 
    }
  }
}
