#!/usr/bin/perl -w
use strict;

##use file from get_insertions_class.pl
my $file = shift;

##use GOSlim from MSU ftp site
my $GO = shift;

my %GO;
my %GO_genes;
open GO, $GO or die "can't open $GO\n";
while ( my $line = <GO> ) {
  next if $line =~ /GOSlim_acc/;
  chomp $line;
  my ( $model, $GOSlim_acc, $GOSlim_name, $aspect ) = split /\t/, $line;
  $GO{$GOSlim_acc}{model}       = $model;
  $GO{$GOSlim_acc}{GOSlim_name} = $GOSlim_name;
  $GO{$GOSlim_acc}{aspect}      = $aspect;
  my ($gene) = $model =~ /^(\S+)\.\d+/;
  ${ $GO_genes{$gene} }{$GOSlim_acc} = 1;
}
print "insert_pos\tgeneID\tinsert_feature\tinsert_type\tGO\n";

open INFILE, "$file" or die;
while ( my $line = <INFILE> ) {
  next if $line eq "\n";
  next if $line =~ /source/;
  chomp $line;
  my (
    $source,   $ref,      $start,          $insert_type,
    $spanners, $flankers, $insert_feature, $f_start,
    $f_end,    $f_strand, $f_name,         $f_note
  ) = split /\t/, $line;
  my $gene = $f_name;
  if ( $insert_feature eq 'intergenic' ) {
    my $dfs = $start - $f_start;
    my $dfe = $f_end = $start;
    ##  LOC_Os10g41820:22503433..22510131(-)/LOC_Os10g41829:22511022..22513298(+)
    my ($g2l,$g2l_strand, $g2r, $g2r_strand) = $f_strand =~ /^(LOC_.+):.+\((.)\)\/(LOC_.+):.+\((.)\)/;
    if ( $dfs < $dfe and $dfs <= 500 ) {
      $gene = $g2l;
      if ($g2l_strand eq '+'){
        $insert_feature = '3prime';
      }else {
        $insert_feature = 'promoter';
      }
      my $go = get_go($gene);
      print "$ref:$start\t$gene\t$insert_feature\t$insert_type\t$go\n";
    }
    elsif ( $dfe < $dfs and $dfe <= 500 ) {
      $gene = $g2r;
      if ($g2r_strand eq '-'){
        $insert_feature = '3prime';
      }else {
        $insert_feature = 'promoter';
      }
      my $go = get_go($gene);
      print "$ref:$start\t$gene\t$insert_feature\t$insert_type\t$go\n";
    }
  }
  elsif ( $insert_feature !~ /mRNA/ ) {
    my $go = get_go($gene);
    print "$ref:$start\t$gene\t$insert_feature\t$insert_type\t$go\n";
  }

}

sub get_go {
  my $gene = shift;
  my @all_GOs;
  foreach my $go_acc ( keys %{ $GO_genes{$gene} } ) {
    my $GO_name   = $GO{$go_acc}{GOSlim_name};
    my $GO_aspect = $GO{$go_acc}{aspect};
    push @all_GOs, "$go_acc, $GO_aspect, $GO_name";
  }
  return join( ";", @all_GOs );
}