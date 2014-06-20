#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $methyl_file = shift;
my $fasta = shift;
my $check_seq = defined $fasta ? 1 : 0;
my %methy;

open FASTA, $fasta or  die "Can't open $fasta";
my %seq;
my $id;
while ( my $line = <FASTA> ) {
  chomp $line;
  if ($line =~ /^>(\S+)/){
    $id = $1;
  }elsif($line !~ /^\s*$/) {
    push @{$seq{ucfirst($id)}} ,  (split ('', $line));
  }
}
open IN, $methyl_file or die "Can't open $methyl_file";
##CHH_RIL16.bismarkextct.out
##CHG_RIL16.bismarkextract.out
my ($type,$strain) = $methyl_file =~ /(C..)\_(.+)\.bismarkext/;
while ( my $line = <IN> ) {
  chomp $line;
  next if $line =~ /Bismark methylation extractor/;
  #DGGXHXP1:377:C2K44ACXX:8:1101:2017:2215_1:N:0:CCGTCC	+	Chr5	11943414	Z
  my (  $read, $strand, $ref, $pos, $methy ) = split /\t/, $line;

  ##A119	CpG	DGGXHXP1:377:C2K44ACXX:8:1101:2017:2215_1:N:0:CCGTCC	+	Chr5	11943414	Z
  ##my ( $strain, $type, $read, $strand, $ref, $pos, $methy ) = split /\t/, $line;

  $methy{$strain}{ucfirst($ref)}{$type}{$pos}{$methy}++;

}
if (!-d "split_methyl"){
mkdir "split_methyl";
}
foreach my $strain ( keys %methy ) {
  foreach my $ref  ( keys %{ $methy{$strain} } ) {
    foreach my $type ( keys %{ $methy{$strain}{$ref} } ) {
      open OUT, ">", "split_methyl/${strain}_${ref}_${type}.split.txt" or die "can't write to ${strain}_${ref}_${type}.split.txt";
      foreach my $pos ( sort {$a<=>$b} keys %{ $methy{$strain}{$ref}{$type} } ) {
        my ( $methyl_yes, $methyl_no ) = ( 0, 0 );
        foreach my $code ( keys %{ $methy{$strain}{$ref}{$type}{$pos} } ) {
          my $count = $methy{$strain}{$ref}{$type}{$pos}{$code};
          if ( $code eq uc $code ) {
            $methyl_yes = $count;
          } else {
            $methyl_no = $count;
          }
        }
        my $nt = uc(${$seq{$ref}}[$pos -1]);
        print OUT join(
                           "\t",
                           (
                             $nt, $pos, $methyl_yes, $methyl_no
                           )
          ),
          "\n";# if uc(${$seq{$ref}}[$pos]) eq 'C';
      }
    }
  }
}
