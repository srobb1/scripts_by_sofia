#!/usr/bin/perl 
use warnings;
use strict;

my $tab = shift;
##CHROM  POS     REF     HEG4_0  HEG4_1  HEG4_2  HEG4_2-5kb
#Chr1    31071   A       G/G     G/G     G/G     G/G
#Chr1    31478   C       T/T     T/T     T/T     T/T
#Chr1    31843   G       A/A     A/A     A/A     A/A
#Chr1    33667   A       G/G     G/G     G/G     G/G
#Chr1    34057   C       T/T     T/T     T/T     T/T

open TAB, $tab or die "Can't open $tab";

##CHROM  POS     REF     EG4_2   HEG4_2  NB
#Chr1    1117    A       A/A     A/C    A/A
chomp (my $header = <TAB>);
my @header = split /\t/ , $header;
shift @header;
shift @header;
shift @header;
my @strains = @header;
my $strain_count = @strains;
print $header , "\n";

while (my $line = <TAB>){
  chomp $line;
  my %SNPs;
  my ($chr , $pos , $ref_nt , @snps) = split /\t/ , $line;
  for (my $i=0 ; $i < $strain_count; $i++){
    my $SNP = $snps[$i];
    next if $SNP eq './.';
    $SNPs{$SNP}++;
  }
  foreach my $SNP (sort { $SNPs{$b} <=> $SNPs{$a}} keys %SNPs ){
    my $count = $SNPs{$SNP};
    my ($a1,$a2) = split /\// ,$SNP;
    next if $a1 eq $a2 and $a1 eq $ref_nt;
    if ($count >= 3){
      print $line,"\n";
    } 
    last;
  }
}
