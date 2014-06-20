#!/usr/bin/perl 
use warnings;
use strict;
use Data::Dumper;


my $inserts = "inserts.near.genes.txt";
my $methy   = "percent_methylation_per_gene.txt";

my %combo;

open INSERTS, $inserts;

#5081	mping	1	Chr1:1041521..1041523	A119_2 Non-reference	LOC_Os01g02890	1	Chr1:1046604..1053166	phosphatidylserine synthase, putative, expressed
#121	mping	1	Chr1:2285174..2285176	A119_2 Non-reference	LOC_Os01g04930	-1	Chr1:2281534..2285053	MYB family transcription factor, putative, expressed
while ( my $line = <INSERTS> ) {
  chomp $line;
  my (
       $dist, $te,       $te_strand, $te_loc, $strain, $insert_type,
       $gene, $g_strand, $g_loc,     $gene_des
  ) = split /\t/, $line;
#  my ( $strain, $insert_type ) = split /\s+/, $strain_type;
  $strain =~ s/(.+)_\d/$1/;
  if ( exists $combo{$gene}{$te_loc}{te_dist}
       and $dist != $combo{$gene}{$te_loc}{te_dist} )
  {
    $dist = 0;
  }

  #  $combo{$gene}{$te_loc}{te_line}=$line;
  $combo{$gene}{$te_loc}{te_strains}{$strain}++;
  $combo{$gene}{$te_loc}{te_dist}   = $dist;
  $combo{$gene}{$te_loc}{te_strand} = $te_strand;

  #  $combo{$gene}{$te_loc}{te_loc}=$te_loc;
  $combo{$gene}{$te_loc}{te_type} = $insert_type;
  $combo{$gene}{$te_loc}{tes}{$te}++;
}
open METHY, $methy;

#gene	ref:start..end	strand	strain	%CpG	%CHH	%CHG	note
#LOC_Os01g04960	Chr1:2299737..2305323	1	A119	92.7	8.0	76.5	transposon protein, putative, unclassified, expressed
      print join( "\t",
                  "TE_in_this_strain", "dist_gene2TE",  "strain",     "gene",
                  "g.loc",     "g.strand",    "g.note",       "g%CpG",
                  "g%CHH",     "g%CHG",       "TE_found_in_Strains", "TEs",
                  "TE_loc",  "TE_strand", "TE_type" ),
        "\n";



while ( my $line = <METHY> ) {
  chomp $line;
  my ( $gene, $loc, $strand, $strain, $CpG, $CHH, $CHG, $note ) = split /\t/,
    $line;
  if ( exists $combo{$gene} ) {
    foreach my $te_loc ( sort keys %{ $combo{$gene} } ) {
#print Dumper $combo{$gene}{$te_loc};
      my ( %te_strains, %tes, $distance, $te_strand, $te_type );
      foreach my $s ( keys %{ $combo{$gene}{$te_loc}{te_strains} } ) {
        $te_strains{$s}++;
      }
      foreach my $t ( keys %{ $combo{$gene}{$te_loc}{tes} } ) {
        $tes{$t}++;
      }
      $distance  = $combo{$gene}{$te_loc}{te_dist};
      $te_strand = $combo{$gene}{$te_loc}{te_strand};
      $te_type   = $combo{$gene}{$te_loc}{te_type};

      my $te_strains = join( '/', sort keys %te_strains );
      my $tes        = join( '/', sort keys %tes );
      my $present    = 'NO';
      if ( exists $te_strains{$strain} ) {
        $present = 'YES';
      }
      print join( "\t",
                  $present, $distance,  $strain,     $gene,
                  $loc,     $strand,    $note,       $CpG,
                  $CHH,     $CHG,       $te_strains, $tes,
                  $te_loc,  $te_strand, $te_type ),
        "\n";
    }
  } else {
    print join( "\t",
                'no_TE', 'no_TE', $strain, $gene,   $loc,
                $strand, $note,   $CpG,    $CHH,    $CHG,
                'no_TE', 'no_TE', 'no_TE', 'no_TE', 'no_TE' ),
      "\n";
  }
}
