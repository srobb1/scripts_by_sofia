#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my @args = @ARGV;
if (!defined @ARGV ){
  die "script.pl file id pingCount file id pingCount file id pingCount\n";
}
my %args;
for (my $i =0; $i < @args ; $i+=3){
  $args{$args[$i]}{id}=$args[$i+1];
  $args{$args[$i]}{ping}=$args[$i+2];
}

my %pos;
foreach my $file (keys %args){
  open IN, $file;
  my $id = $args{$file}{id};
  my $ping = $args{$file}{ping};
  while (my $line = <IN>){
    chomp $line;
    my $diff=0;
    my ($pos, $A, $T, $G, $C, $N, $total)= split /\t/, $line; 
      foreach my $nt_count ($A, $T, $G, $C ){
        my ($nt , $count) = split /\:/ , $nt_count;
        #if ($count > ($total/$ping * .2)){
        if ($count > ($total/$ping * 0.5)){
          $diff++;
        }
        if ($nt eq 'A'){
          $A =  $count;
        }if ($nt eq 'T'){ 
          $T =  $count;
        }if ($nt eq 'G'){ 
          $G =  $count;
        }if ($nt eq 'C'){ 
          $C =  $count;
        }
      }
      if ($diff > 1){
        $pos{$pos}{$id}{A}=$A;
        $pos{$pos}{$id}{T}=$T;
        $pos{$pos}{$id}{G}=$G;
        $pos{$pos}{$id}{C}=$C;
        $pos{$pos}{$id}{total}=$total;
      }
  }
}

foreach my $pos (sort { (split /\(/, $a)[0] <=> (split /\(/ ,$b)[0]} keys %pos){
  my $line = "$pos ::";
  my $diffcount = 0;
  my $HEG4_count= 0;
  foreach my $id (sort keys %{$pos{$pos}}){
    $HEG4_count++ if $id =~ /HEG4/;
    my $total = $pos{$pos}{$id}{total};
    $line .= "\n\t$id($total)\t";
      $diffcount++;
    foreach my $nt (sort keys %{$pos{$pos}{$id}}){
      next if $nt eq 'total';
      my $count = $pos{$pos}{$id}{$nt};
      $line .=  "$nt:$count\t";
    }
  }
  $line =~s/\t$/\n/;
  $line =~s/\s::/:($diffcount)/;
  next if $HEG4_count == 1;
  next if $line =~ /HEG4_0/ and $diffcount eq 1;
  print $line;
}
