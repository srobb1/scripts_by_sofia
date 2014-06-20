#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;
my $seqIO_obj = Bio::SeqIO->new(-file=>$file, -format=>"fasta");
open STUDENT ,">convertedStudent.fa" or die "Can't open convertedStudent.fa";
open INSTRUCTOR ,">convertedInstructor.fa" or die "Can't open convertedInstrucot.fa";
while (my $seq_obj = $seqIO_obj->next_seq){
  my $id = $seq_obj->id;
  my $seq = lc($seq_obj->seq);
  my @seq = split '',$seq;
  my @Cs;
  for (my $i=0 ; $i<@seq ; $i++){
    my $nt = $seq[$i];
    if ($nt eq 'c'){
      push @Cs,$i;
    } 
  }
  my ($min,$max) = (10,60);
  my $percent_whole = int(rand(($max - $min) + 1)) + $min;
  my $frac = rand;
  $percent_whole += $frac;
  my $percent = $percent_whole/100;
  my $C_count = scalar @Cs;
  my $convertedCs = int((scalar @Cs) * $percent);
  my %alreadyConverted;
  for (my $i = 0 ; $i < $convertedCs ; $i++){
    my $rand =  int(rand($C_count));
    my $pos = $Cs[$rand];
    #my $sl = length $seq;
    #print "$id($sl): $percent_whole\% is $convertedCs $rand $pos\n";
    if (!exists $alreadyConverted{$pos}){
      $seq[$pos] = 'T';
      $alreadyConverted{$pos}++;
    }else{
      $i--;
    }

  }
  my $convertedSeq = join '', @seq;
  my $C_orig = $seq =~ tr/c/c/; 
  my $C_convert = $convertedSeq =~ tr/T/T/; 
  my $per_convert = $C_convert/$C_orig;
  print INSTRUCTOR ">$id $percent_whole\%\nORIG:$seq\nCONV:$convertedSeq\n";
  $convertedSeq =~ tr/T/t/;
  print STUDENT ">$id\n$convertedSeq\n";
}
