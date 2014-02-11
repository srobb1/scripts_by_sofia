#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my @gffs = @ARGV;

my %inserts;
foreach my $gff (@gffs){
  open IN, $gff or die "Can't open $gff for reading \n";
  while (my $line = <IN>){
    next if $line =~ /^#/;
    my ($ref,$strain,$type,$start,$end,$score,$strand,$phase,$nine)=split /\t/ , $line;
    my ($TE) = $nine =~ /ID=$strain\.(.+)\.te_insertion/;
    my ($left,$right) = $nine =~ /left_flanking_read_count=(\d+);right_flanking_read_count=(\d+)/;
    my $st = $strand !~ /[+-]/ ? 'UNK' : $strand;
    if (defined $left and defined $right){
      $inserts{$ref}{$start}{$end}{$strain}{$TE}="($st)$left|$right";
    }else{
      $inserts{$ref}{$start}{$end}{$strain}{$TE}="($st)";
    }
  }
 }
#print Dumper \%inserts;
foreach my $ref (keys %inserts){
  foreach my $start (keys %{$inserts{$ref}}){
    foreach my $end (keys %{$inserts{$ref}{$start}}){
      my @inserts;
      foreach my $strain(keys %{$inserts{$ref}{$start}{$end}}){
        my @TEs;
        foreach my $TE(keys %{$inserts{$ref}{$start}{$end}{$strain}}){
          my $info= $inserts{$ref}{$start}{$end}{$strain}{$TE};
          push @TEs, "$TE$info"; 
        }
        push @inserts , "$strain:". join (':',@TEs);
      }
      my $inserts = "inserts=".join (',',@inserts);
      my $ID = "ID=te_insertion_site.$ref.$start.$end";
      print join("\t",$ref,'.','transposable_element_insertion_site',$start,$end,'.','.','.',"$ID;$inserts"),"\n";
    }
  }
}

