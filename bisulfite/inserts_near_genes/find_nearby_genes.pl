#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;

my $sqlite = shift;
die "Please provide a SQLite datafile of a seqfeature db" if !defined $sqlite;

open BED , ">inserts_and_genes.bed" or die "Can't Open gene_with_upstream_retro.txt for writting\n";

#############################################################################
## make a Bio::DB::SeqFeature::Store object (contains info about  organism in
## SQLite database and all the methods needed to interact with the database
#############################################################################

#Open the sequence database
my $db_obj = Bio::DB::SeqFeature::Store->new(
  -adaptor => 'DBI::SQLite',
  -dsn     => $sqlite
);

## get_feature_by_type can get any type from your original GFF, anything in Col3
my @features_type_TE = $db_obj->get_features_by_type('transposable_element_insertion_site');
#my @features_type_TE = $db_obj->get_features_by_location('Chr9',10863118, 10863120);
foreach my $feature (@features_type_TE) {
  my $id  = $feature->load_id;
  my @id = split /\./ , $id;
  my $TE_name = $id[1];
  next unless $feature->type =~ /transposable_element_insertion_site/;
  my $TE_start = $feature->start;
  my $name = $feature->name;
  my $TE_end   = $feature->end;
  my $ref        = $feature->ref;
  my $strand = $feature->strand;
  my $source = $feature->source;
  my %attr = $feature->attributes;
  my $note = ${$attr{Note}}[0];
  if ($id =~ /\.ping\.te_insertion_site/){
    #$note = "$name $source"; 
    $note = "$source\t$name"; 
  }elsif($id =~ /NB\..+\.te_insertion_site/){ 
    #$note = "$TE_name $source";
    $note = "$source\t$note";
  }else{
    if (!defined $note){
     #print "$id\n";
     #print Dumper \%attr;
    }else{
      $note = $id[0] . "\t$note";
    }
  }
  #my $s = $strand > 0 ? '+' : '-';
  my $TE_info = "$TE_name\t$strand\t$ref:$TE_start..$TE_end\t$note";
  my @dists_from_TE = (10000);
  ## find each TE, look upstream 1000, then 2000, then 5000. stop at first gene feature found
  foreach my $dist(@dists_from_TE){
    ## now look upstream of TEs on the plus and minus strand
    my @coords = sort {$a <=> $b }  ( $TE_start, $TE_start-$dist );
    my @features_left  =  $db_obj->get_features_by_location($ref,$coords[0], $coords[1]);
    @coords = sort {$a <=> $b }  ( $TE_end, $TE_end+$dist );
    my @features_right =  $db_obj->get_features_by_location($ref,$coords[0], $coords[1]);

    ## for genes 2 the left
    foreach my $f ( sort { $b->end <=> $a->end } @features_left){
      my $type = $f->type;
      next if $type =~ /transposable_element_insertion_site/;
      next unless $type =~ /gene/;
      my $gene_start = $f->start;
      my $gene_end   = $f->end;
      my $strand = $f->strand;
      my %attr= $f->attributes;
      my $note = ${$attr{'Note'}}[0];
      my $gene_name = $f->name;
      my $gene_info = "$gene_name\t$strand\t$ref:$gene_start..$gene_end\t$note";
      my $distance = $TE_start - $gene_end;
      print "$distance\t$TE_info\t$gene_info\n";
      last;
    }
    ## do it again for genes 2 the right
    foreach my $f ( sort { $a->end <=> $b->end } @features_right){
      my $type = $f->type;
      next if $type =~ /transposable_element_insertion_site/;
      next unless $type =~ /gene/;
      my $gene_start = $f->start;
      my $gene_end   = $f->end;
      my $strand = $f->strand;
      my %attr= $f->attributes;
      my $note = ${$attr{'Note'}}[0];
      my $gene_name = $f->name;
      my $gene_info = "$gene_name\t$strand\t$ref:$gene_start..$gene_end\t$note";
      my $distance = $gene_start - $TE_end;
      print "$distance\t$TE_info\t$gene_info\n";
      last;
    }

  }
}
