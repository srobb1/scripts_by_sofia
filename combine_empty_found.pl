#!/usr/bin/perl -w
use strict;
use Data::Dumper;

## perl combine_empty_found.pl 08202012.parse_blat_empty.summary.NAM.txt  08072012.NAM.parsed.out  dups_te.list

#==> 08202012.parse_blat_empty.summary.NAM.txt <==
#library te      confirmedEnds   read_count
#P39     78|ZM_PIF/Harbinger_22|chr10:85500875..85501053 empty_site      1
#
#==> 08072012.NAM.parsed.out <==
#library te      confirmedEnds   five_count      three_count
#CML333  1|ZM_hAT_noncoding_170|chr10:1134757..1134821   one-end 3       0
#
#==> dups_te.list <==
#1091|ZM_PIF/Harbinger_26|chr2:234432781..234433045|five_prime|chr2:234432756..234432805
#1097|ZM_hAT_noncoding_42|chr2:235921482..235921658|five_prime|chr2:235921457..235921506

my $empty = shift;
my $found = shift;
my $dups = shift;

my %dups;
open DUPS, $dups or die "Can't open $dups";
while (my $line = <DUPS>){
  #2786|ZM_Tourist_34|methyl|chr5:143669687..143669802|complete|chr5:143669687..143669802|five_prime|chr5:143669662..143669711
  chomp $line;
  my ($num_id,$te_id,$three,$methyl_range,$five,$range,$desc,$desc_range) = split /\|/ , $line;
  my $id = "$num_id|$te_id|$range";
  $dups{$id} = $line;
}

open EMPTY , $empty or die "Can't open $empty";
my %empty;
while (my $line = <EMPTY>){
  next if $line =~ /library/;
  chomp $line;
  #P39     78|ZM_PIF/Harbinger_22|chr10:85500875..85501053 empty_site      1
  my ($NAM , $id, $desc,$count) = split /\s+/ , $line;
  $empty{$NAM}{$id} = $count;
}
#print Dumper \%empty;


close EMPTY;
open FOUND, $found or die "Can't open $found";
print join ("\t", qw(NAM TE_id 5'readCount 3'readCount emptySiteCount)) , "\n"; 

while (my $line = <FOUND>){
  #CML333  1|ZM_hAT_noncoding_170|chr10:1134757..1134821   one-end 3       0
  next if $line =~ /library/;
  chomp $line;
  my ($NAM, $id, $ends, $one, $two) = split /\s+/ , $line;
  my $empty = 0;
  if (exists  $empty{$NAM}{$id} ){
    $empty = $empty{$NAM}{$id};
    delete $empty{$NAM}{$id};
  }
  my $dup = exists $dups{$id} ? $dups{$id} : 0;
  if ($dup =~ /five/) {
    $one = $one.".dup";
  } elsif ($dup =~ /three/){
    $two = $two.".dup";
  }
  print join ("\t", $NAM, $id, $one, $two, $empty) , "\n"; 
}
foreach my $NAM (sort keys %empty){
  foreach my $id (sort keys %{$empty{$NAM}}){
   my $empty = $empty{$NAM}{$id};
   print join ("\t", ($NAM,$id,0,0,$empty)) , "\n";
  }

}
