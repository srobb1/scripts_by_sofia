#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Bio::DB::Fasta;
my $ortho_table = shift; #mclGroups.I15.table
my $good_proteins_fasta = shift; #goodProteins.fasta

my $db      = Bio::DB::Fasta->new($good_proteins_fasta);

open IN, $ortho_table or die "Can't open $ortho_table\n";
my %orthologs;

my @ids = `grep '>' $good_proteins_fasta`;

while (my $line=<IN>){
#OG7389: MorElo_Prot|882244 MorVer_Prot|MVEG_05998
  chomp $line;
  my ($ortho, $proteins) = split /\:/ , $line;
  my @proteins = split ' ' ,$proteins;
  foreach my $protein(@proteins){
    $orthologs{$protein}=$ortho;
  }
}
#print Dumper \%orthologs;
foreach my $line (@ids){
    chomp $line;
    my $id = substr $line , 1;
    
    if (! exists $orthologs{$id}){
      #print $id,"\n";
      my $seq_obj = $db->get_Seq_by_id($id);
      my $seq = $seq_obj->seq;
      print ">$id\n$seq\n";
    }
}

__END__
foreach my $group (keys %orthologs){
  my $count = scalar keys %{$orthologs{$group}};
  if ($count == 2){
   my ($spec) =  keys %{$orthologs{$group}};
   print "$group\t$spec\n";
  }
}
