#!/usr/bin/perl -w
use strict;
use Data::Dumper;

## #CHROM  POS     REF     HEG4_2  NB      RIL39
## Chr2    17321   A       G/G     A/A     A/A
## Chr2    60863   T       G/G     T/T     T/T
## Chr2    81426   C       T/T     C/C     C/C
## Chr2    114497  C       T/T     C/C     C/C
## Chr2    118723  C       G/G     C/C     C/C
## Chr2    122239  A       A/G     A/A     A/A
## Chr2    123298  G       A/A     G/G     G/G
## Chr2    133676  C       G/G     C/C     C/C
## Chr2    390336  T       C/C     T/T     T/T
## Chr2    427338  C       T/T     C/C     C/C
## Chr2    441152  T       T/C     T/T     T/C
## Chr2    474264  G       G/A     G/G     G/A

my $file = shift;
open INFILE, $file or die "cant open $file\n";
my $header = <INFILE>;
chomp $header;
my ($h0 , $h1, @strains) = split  /\t/ , $header;

my %strains;
for (my $i=0 ; $i < @strains ; $i++){
  $strains{$strains[$i]}=$i;
}

print "chr\tpos\tvalue\n";
while (my $line = <INFILE>){
  my %SNPs;
  chomp $line;
  my ($chr, $pos, @snps) = split "\t" , $line;
  next if $chr eq 'ChrUn' or $chr eq 'ChrSy';
  for (my $i=0 ; $i<@strains ; $i++){
    my $a = $snps[$i];
    my $strain = $strains[$i];
    next if $strain eq 'REF';
    #my @alleles = split '/' , $snp_str;
    #foreach my $a (@alleles){
      $SNPs{$a}{$strain} = 1;
    #}  
  }
  my @alleles = sort keys %SNPs;
  #next if @alleles > 2;
  
  foreach my $a (sort keys %SNPs){
    my @strains = sort keys %{$SNPs{$a}};
    #print "$pos $a @strains\n";
    if (@strains == 3){
       #print "$chr\t$pos\t0\n";
    }elsif(@strains == 2 and $strains[0] eq 'HEG4_2' and $strains[1] eq 'RIL39'){
       #print "$chr\t$pos\t-1\n";
       print "HEG4,RIL39\t$chr.$pos\tsnp\n";
    }elsif(@strains == 2 and $strains[0] eq 'NB' and $strains[1] eq 'RIL39'){
       print "NB,RIL39\t$chr.$pos\tsnp\n";
       #print "$chr\t$pos\t1\n";
    }#else dont print anything
    
  }
}
