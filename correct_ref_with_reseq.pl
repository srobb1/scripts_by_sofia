#!/usr/bin/perl -w
use strict;
use Data::Dumper;
## run 'cat file.vcf | vcf-to-tab > out.tab' before running this script

my $vcf_tab = shift;
open VCFTAB, "$vcf_tab" or die;


##CHROM  POS     REF     EG4_2   HEG4_2  NB
#Chr1    1117    A       A/A     A/C    A/A
chomp (my $header = <VCFTAB>);
my @header = split /\t/ , $header;
shift @header;
shift @header;
shift @header;
my @strains = @header;
my $strain_count = @strains;
$header =~ s/\tNB//;
print $header , "\n";
while (my $line = <VCFTAB>){
  chomp $line;
  my ($chr , $pos , $ref_nt , @snp) = split /\t/ , $line;
  #print "$chr , $pos , $ref_nt , @snp\n";
  my $reseq_index;
  for (my $i=0 ; $i < $strain_count ; $i++){
    if ( $header[$i] eq 'NB'){
      $reseq_index =  $i  ;
    }
  }
  my ($a1 , $a2) = $snp[$reseq_index] =~ /(.)\/(.)/;
  #change ref if diff from reseq NB
  if ( ($a1 eq $a2) and ($a1 ne $ref_nt) and $a1 ne '.'){
    $ref_nt = $a1;
  }
  #what should i do if reseq NB is het at a position
  #throw out any lines that have reseq NB as heter
  elsif ($a1 ne $a2){
    #move to next position, dont print this one    
    next;
  }
  
  print "$chr\t$pos\t$ref_nt";
  for (my $i=0 ; $i < $strain_count ; $i++){
    next if $i == $reseq_index;
    print "\t$snp[$i]";
  }
  print "\n";
}
