#!/usr/bin/perl -w
use strict;

my $vcf = shift;
##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  A119_2  A123_2  EG4_2   HEG4_2  NB
#Chr1    31071   .       A       G       9183.67 PASS    AC=8;AF=0.80;AN=10;BaseQRankSum=1.398;DP=330;Dels=0.00;FS=2.934;HRun=0;HaplotypeScore=0.8109;MQ=59.42;MQ0=0;MQRankSum=1.652;QD=36.44;ReadPosRankSum=-1.126      GT:AD:DP:GQ:PL  1/1:0,42:42:99:1371,99,0        1/1:0,54:54:99:1889,138,0       1/1:0,76:76:99:2959,223,0       1/1:0,80:80:99:2964,223,0       0/0:78,0:78:99:0,166,2303

## GT
##0/0 - the sample is homozygous reference
##0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
##1/1 - the sample is homozygous alternate In the three examples above, NA12878 is observed with the allele combinations T/G, G/G, and C/T respectively.

## DP field describes the total depth of reads that passed the caller's internal quality control metrics (like MAPQ > 17, for example), the 
## AD values (one for each of REF and ALT fields) is the unfiltered count of all reads that carried with them the REF and ALT alleles.

open IN, $vcf or die "Can't open $vcf\n";
my @header;
while (my $line = <IN>){
  my %strains;
  chomp $line;
  my ($CHROM,$POS,$ID,$REF,$ALT,$QUAL,$FILTER,$INFO,$FORMAT,@strains) = split /\t/, $line;
  if ($line =~ /^#CHROM\tPOS/){
    @header = @strains;
    next;
  }elsif($line =~ /^#/){
    next;
  }else{
    for  (my $i = 0 ; $i < @header ;$i++){ 
      my $strain_data = $strains[$i];
      my @format_headers = split /:/ , $FORMAT;
      my @format_data = split /:/ , $strain_data;
      ## 1/1:0,42:42:99:1371,99,0
      if ($strain_data ne './.'){
      for (my $j=0 ; $j < @format_headers ; $j++){
        if ($format_headers[$j] eq 'GT'){
          my $GT = $format_data[$j];
          my $genotype;
          my ($r_freq,$a_freq);
          if ($GT eq '0/0'){
            $genotype = "$REF/$REF";
            $r_freq = 1;
            $a_freq = 0;
          }elsif($GT eq '0/1'){
            $genotype = "$REF/$ALT";
            $r_freq = 0.5;
            $a_freq = 0.5;
          }elsif($GT eq '1/1'){
            $genotype = "$ALT/$ALT"; 
            $r_freq = 0;
            $a_freq = 1;
          }
          $strains{$header[$i]}{GT}=$genotype; 
          $strains{$header[$i]}{ref_freq}=$r_freq; 
          $strains{$header[$i]}{alt_freq}=$a_freq; 
        }elsif ($format_headers[$j] eq 'AD'){
          my ($ref_count,$alt_count) = split /,/ , $format_data[$j];
          $strains{$header[$i]}{ref_count}=$ref_count;
          $strains{$header[$i]}{alt_count}=$alt_count;
        }
      }
      }else{
          $strains{$header[$i]}{ref_count}=0;
          $strains{$header[$i]}{alt_count}=0;
          $strains{$header[$i]}{ref_freq}=0;
          $strains{$header[$i]}{alt_freq}=0;
      }
    }
  }
  my @acounts;
  my @alleles;
  foreach my $strain (keys %strains){
    my $ref_count = $strains{$strain}{ref_count};
    my $alt_count = $strains{$strain}{alt_count};
    my $alt_freq = $strains{$strain}{alt_freq};
    my $ref_freq = $strains{$strain}{ref_freq};
    my $total = $alt_count + $ref_count;
    push @acounts , "$strain:$REF+$ref_freq+$ref_count+$ALT+$alt_freq+$alt_count+$total+GATK";
  }
  my $acount = join ',' , @acounts;
  my $alleles = "$REF/$ALT";
  print join ("\t", $CHROM,"GATK",'snp',$POS,$POS,'.','+','.',"ID=SNP.$CHROM.$POS;Name=SNP.$CHROM.$POS;acount=$acount;alleles=$alleles;refallele=$REF"),"\n";
}
