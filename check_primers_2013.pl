#!/usr/bin/perl -w
use strict;


my $file = shift;
my $db = shift;
my $uniq_only = defined $ARGV[2] ? $ARGV[2] : 0 ;
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
  $primers{$id}{found}=0; 
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
  next if  ($mm > 0);
  next if  ($gapOpenings > 0);
  my $strand = ($qEnd - $qStart > 0) ? 1 : -1 ;
  my $p;
  if ($qStart < $p1_len) {
    ## this is p1 aligned
    next if $alignLen < $p1_len;
    $p = 'p1';
  }else{
    ## ths is p2 aligned
    next if $alignLen < $p2_len;
    $p = 'p2';
  }
  ## mayabe i can use strand and p1 p2 identity to simplify proper primer finding
  $primers{$primer_pair}{blast}{$count}{$subject}{strand} = $strand;
  $primers{$primer_pair}{blast}{$count}{$subject}{alnLen} = $alignLen;
  $primers{$primer_pair}{blast}{$count}{$subject}{perID}  = $perId;
  $primers{$primer_pair}{blast}{$count}{$subject}{qStart} = $qStart;
  $primers{$primer_pair}{blast}{$count}{$subject}{qEnd}   = $qEnd;
  $primers{$primer_pair}{blast}{$count}{$subject}{sStart} = $sStart;
  $primers{$primer_pair}{blast}{$count}{$subject}{sEnd}   = $sEnd;
  $primers{$primer_pair}{blast}{$count}{$subject}{e}      = $e;
  $primers{$primer_pair}{blast}{$count}{$subject}{line}   = $line;
  $primers{$primer_pair}{subject}{$subject}=1;
  $primers{$primer_pair}{found}=1;
  
}
my %warn;
  print "primer_set\tproductSize\thit:start..end\tin_range\tmaps\tprimer1\tprimer2\n";
my %matches;
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
  my ($loc_chr,$loc_start,$loc_end) = (0,0,0);
  if ($primer_set =~ /\b(\w+):(\d+)..(\d+)/){
    ($loc_chr,$loc_start,$loc_end) = ($1, $2, $3);
  }
  foreach my $sub (keys %hits){
    my @sorted_values = sort {$a <=> $b} @{$hits{$sub}};
    my $maps;
    
    if (@sorted_values < 4){
      #2 starts, 2 ends ## just needs to map completely
      # ==========       ===========
      # S        E       S         E
      $maps = "both primers do not map";
    }elsif (@sorted_values == 4){
      #2 starts, 2 ends ## exactly maps 
      $maps = "both primers map uniquely";
    }elsif (@sorted_values > 4){
      #more than only 2 starts and  2 ends  
      $maps = "primers do not map uniquely";
    }
    if ($uniq_only){
      next if @sorted_values != 4; #exactly 2 starts, 2 ends ## needs to map uniquely
    }
    my $smallest = shift @sorted_values;
    my $largest = pop @sorted_values;
    my $p1 = uc $primers{$primer_set}{p1}{seq};
    my $p2 = uc $primers{$primer_set}{p2}{seq};
    my $product_size = $largest-$smallest+1;
    my $in_range = 'N/A';
    if ($sub eq $loc_chr and $loc_start > $smallest and $loc_end < $largest){
      $in_range = 'yes';
    }else{
      $in_range = 'no';
    }
      print "$primer_set\t$product_size\t$sub:$smallest..$largest\t$in_range\t$maps\t$p1\t$p2\n";# if $product_size < 10000 and $product_size > 100;

    #print "$primer_set\t",$largest-$smallest+1,"\t$subjects[0]:$smallest..$largest\t$p1\t$p2\n";
  }
}
foreach my $id (keys %primers){
  my $found = $primers{$id}{found};
  if (!$found){
    print "$id\tNo blast hit found\n";
  }
}
