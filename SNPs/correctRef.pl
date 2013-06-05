#!/usr/bin/perl 
use strict;
use warnings;
use Bio::SeqIO;

### the $REF contains the value the individual you would like to use to replace the
### original reference nucleotide

my $VCF_tab = shift;
my $genome = shift;
my $REF = 'NB';
my $seqIO_obj = Bio::SeqIO->new (-file=>$genome , -format=>'fasta');
open VCFTAB , $VCF_tab or die "Can't Open $VCF_tab\n";

chomp ( my $header = <VCFTAB> );
##CHROM  POS     REF     HEG4_0  HEG4_1  HEG4_2  HEG4_2-5kb      NB
#Chr1    31071   A       G/G     G/G     G/G     G/G     A/A
my ($one, $two, $three, @strains) =  split /\t/ , $header;
my $REF_i = '';
for ( my $i=0 ; $i<@strains ; $i++){
  if ($strains[$i] eq $REF){
    $REF_i = $i;
  }
}
my %correction;
while (my $line = <VCFTAB>){
  chomp $line;
  my ($seq, $position, $ref_nt, @SNPs) = split /\t/ , $line;
  my $REF_snp = $SNPs[$REF_i];
  my ($a1,$a2) = $REF_snp =~ /(.)\/(.)/;
  next if $a1 eq '.';
  next if $a1 ne $a2;
  if ($ref_nt ne $a1){
   $correction{$seq}{$position}=$a1;
  }
}
my $out_seqIO = Bio::SeqIO->new(-format=>"fasta");
while ( my $seq_obj = $seqIO_obj->next_seq){
  my $id = $seq_obj->id;
  my $seq = $seq_obj->seq;
  my @seq = split '' , $seq;
  unshift @seq, '';
  foreach my $position (sort {$a <=> $b} keys %{$correction{$id}}){
    my $nt = $correction{$id}{$position};
    $seq[$position] = $nt;
  }
  shift @seq ; # remove the 0th element
  my $newseq = join('',@seq);
  $newseq =~ s/(.{80})/$1\n/g; 
  print ">$id\n$newseq\n";
}

