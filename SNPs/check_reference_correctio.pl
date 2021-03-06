#!/usr/bin/perl
use strict;
use warnings;
use Bio::DB::Fasta;

my $org = shift;
my $new = shift;
my $vcf_tab = shift;


my $genome = shift;
my $REF = 'NB';
open VCFTAB , $vcf_tab or die "Can't Open $vcf_tab\n";

my $new_DBobj = Bio::DB::Fasta->new ($new);
my $org_DBobj = Bio::DB::Fasta->new ($org);

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
   $correction{$seq}{$position}{new}=$a1;
   $correction{$seq}{$position}{old}=$ref_nt;
  }
}

foreach my $id (keys %correction){
  my $seq_new = $new_DBobj->seq($id);
  my $seq_org = $org_DBobj->seq($id);
  my @seq_new = split '' , $seq_new;
  my @seq_org = split '' , $seq_org;
  unshift @seq_new , '';
  unshift @seq_org , '';
  foreach my $position ( sort {$a <=> $b} keys %{$correction{$id}}){
    my $new_nt = $correction{$id}{$position}{new};
    my $org_nt = $correction{$id}{$position}{old};
    print "$id:$position\t$org_nt\t$seq_org[$position]\t$new_nt\t$seq_new[$position]\n";

  }  
}
