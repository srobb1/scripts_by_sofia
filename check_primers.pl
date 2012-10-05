#!/usr/bin/perl -w
use strict;


my $file = shift;
my $db = shift;

my %primers;

if (!defined $file){
  die "please provide a fasta file of primers in the following format:
primersetName PRIMER_foward PRIMER_rev
primersetName PRIMER_foward PRIMER_rev
primersetName PRIMER_foward PRIMER_rev
";
}

if (!defined $db){
  die "please provide a fasta to use as a db

";
}
open IN , $file or die "Can't open $file $!";
open OUT , ">$file.fa" or die "Can't open $file.fa for writting $!";

while (my $line = <IN>){
  my ($id , $primer1,$primer2 ) = split /\s+/ , $line;
  $primers{$id}{p1}{seq}=$primer1;
  $primers{$id}{p2}{seq}=$primer2;
  $primers{$id}{p1}{len}=length ($primer1);
  $primers{$id}{p2}{len}=length ($primer2);
  my $seq_for_blast = $primer1 . "NNNNNNNNNN" . $primer2;
  
 print OUT ">$id\n$seq_for_blast\n";

}


`blastall -F F -p blastn -i $file.fa -d $db -o $file.blastout -m8` if !-e "$file.blastout";


my $count;
open BLASTOUT, "$file.blastout" or die "Can't open $file.blastout";
#1	Chr11	100.00	21	0	0	45	65	21964649	21964629	0.004	42.1
#1	Chr11	100.00	16	0	0	1	16	21963880	21963895	3.6	32.2
while ( my $line = <BLASTOUT> ) {
  $count++;
  chomp $line;
  my (
    $primer_pair, $subject,     $perId,  $alignLen,
    $mm,          $gapOpenings, $qStart, $qEnd,
    $sStart,      $sEnd,        $e,      $score
  ) = split /\t/, $line;
  my $p1_len = $primers{$primer_pair}{p1}{len};
  my $p2_len = $primers{$primer_pair}{p2}{len};
  next if  (($alignLen != $p1_len) and ($alignLen != $p2_len));
  next if  ($perId != 100);
  #$primers{$primer_pair}{blast}{$count}{sub}    = $subject;
  $primers{$primer_pair}{blast}{$count}{$subject}{alnLen} = $alignLen;
  $primers{$primer_pair}{blast}{$count}{$subject}{mm}     = $mm;
  $primers{$primer_pair}{blast}{$count}{$subject}{perID}  = $perId;
  $primers{$primer_pair}{blast}{$count}{$subject}{qStart} = $qStart;
  $primers{$primer_pair}{blast}{$count}{$subject}{qEnd}   = $qEnd;
  $primers{$primer_pair}{blast}{$count}{$subject}{sStart} = $sStart;
  $primers{$primer_pair}{blast}{$count}{$subject}{sEnd}   = $sEnd;
  $primers{$primer_pair}{blast}{$count}{$subject}{e}      = $e;
  $primers{$primer_pair}{blast}{$count}{$subject}{line}   = $line;
  $primers{$primer_pair}{subject}{$subject}=1;
  
}
my %warn;
  print "primer_set\tproductSize\thit:start..end\tin_range\tprimer1\tprimer2\n";
foreach my $primer_set ( keys %primers ) {
  #my @subjects = keys %{$primers{$primer_set}{subject}};
  #next if @subjects > 1;
  my %hits;
  my @values;
  my @lines;
  foreach my $count (keys %{$primers{$primer_set}{blast}}){
    foreach my $sub (keys %{$primers{$primer_set}{blast}{$count}}){
      push @{ $hits{$sub}} , $primers{$primer_set}{blast}{$count}{$sub}{sStart} , $primers{$primer_set}{blast}{$count}{$sub}{sEnd};
    }
  }
  my $loc = 0;
  if ($primer_set =~ /\D+(\d{3,})/){
    $loc = $1;
  }
  foreach my $sub (keys %hits){
    my @sorted_values = sort {$a <=> $b} @{$hits{$sub}};
    next if @sorted_values < 4; #2 starts, 2 ends
    my $smallest = shift @sorted_values;
    my $largest = pop @sorted_values;
    my $p1 = uc $primers{$primer_set}{p1}{seq};
    my $p2 = uc $primers{$primer_set}{p2}{seq};
    my $product_size = $largest-$smallest+1;
    my $in_range = 'N/A';
    if (defined $loc and $loc > $smallest and $loc < $largest){
      $in_range = 'yes';
    }else{
      $in_range = 'no';
    }
    print "$primer_set\t$product_size\t$sub:$smallest..$largest\t$in_range\t$p1\t$p2\n" if $product_size < 10000 and $product_size > 100;
    #print "$primer_set\t",$largest-$smallest+1,"\t$subjects[0]:$smallest..$largest\t$p1\t$p2\n";
  }
}
