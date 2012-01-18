#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# Open the sequence database
my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::SQLite',
                                           -dsn     => '/home_stajichlab/robb/rice/database/MSUr7.mping.insertions.sqlite');
my %ref;
my $ref_total = 0;
my %transposons;
my @features_type = $db->get_features_by_type('chromosome');
foreach my $feature (@features_type) {
  my $ref = $feature->ref;
  next if $ref =~/Un/ or $ref =~ /Sy/;
  my $len = $feature->length;
  $ref{$ref}= $len ;
  $ref_total+=$len;
}
print "ref: $ref_total\n";
my %feat;
my %strand;
foreach my $feature ($db->get_features_by_type('gene')) {
  my $ref = $feature->ref; 
  if (!exists $feat{intergenic}{$ref}){
    for (my $i = 0; $i < $ref{$ref} + 1 ; $i++){
      ## set every bp to 1
      ${$feat{intergenic}{$ref}}[$i] = 1;
      ## to keep track of gene strand
      ${$strand{$ref}}[$i] = 0;
      ${$feat{intergenic}{$ref}}[0] = 0;
    }
  }

### do not count anything that belongs to a gene that is annotated as a transposon
  my %attr = $feature->attributes;
  my $notes_arry_ref = $attr{'Note'};
  my $note = ${$notes_arry_ref}[0];
  my $name = $feature->seq_id;
  next if $note =~ /transposon/;
## end transposon check

  my $f_start = $feature->start;
  my $f_end = $feature->end;
  ## keep track of gene strand for intergenic features
  my $f_strand = $feature->strand;
  for (my $i= $f_start ; $i < ($f_end +1) ; $i++){
    ${$strand{$ref}}[$i] = $f_strand; 
  }
  ## end intergenic feature strand

  ## where a gene exists change intergenic index to 0
  for (my $i= $f_start ; $i < ($f_end +1) ; $i++){
    ${$feat{intergenic}{$ref}}[$i] = 0; 
  }
}  

  my $feat_total = 0;
  foreach my $ref (sort keys %{$feat{intergenic}}){
    my $feat_sum = sum @{$feat{intergenic}{$ref}}; 
    $feat_total += $feat_sum;
  }
  print "intergenic: $feat_total\n";
my %intergenic;
foreach my $ref (sort keys %ref){
my $last_value = 0;
my $start;
my $end;
my $gene_2_right_strand;
my $gene_2_left_strand;

  for (my $i=0; $i < $ref{$ref} + 1 ; $i++){
    my $value = ${$feat{intergenic}{$ref}}[$i];
      if ($last_value == 0 and $value == 1){
        ##first position of intergenic
        $start = $i;
        $gene_2_left_strand = $strand{$ref}[$i-1]
       }elsif($last_value == 1 and $value == 0) {
        ##last position of intergenic
        $end = $i-1;
        $gene_2_right_strand = $strand{$ref}[$i];
        $intergenic{$ref}{$start}{end} = $end;
        $intergenic{$ref}{$start}{left_strand} = $gene_2_left_strand;
        $intergenic{$ref}{$start}{right_strand} = $gene_2_right_strand;
        ## reset variables
        $gene_2_right_strand = '';
        $gene_2_left_strand = '';
	$start = '';
        $end = '';
      }
      $last_value = $value;
  }
}
open GFF, ">intergenic.gff" or die "Can't open intergenic.gff for writing\n";
foreach my $ref (sort keys %intergenic){
  foreach my $start (sort {$a <=> $b} keys %{$intergenic{$ref}}){
    my $end = $intergenic{$ref}{$start}{end};
    my $gene_2_left_strand = $intergenic{$ref}{$start}{left_strand};
    my $gene_2_right_strand = $intergenic{$ref}{$start}{right_strand};
    my $len = $end - $start + 1;
    print GFF "$ref\t.\tintergenic\t$start\t$end\t.\t.\t.\tID=$ref.intergenic.$start.$end;gene_2_left_strand=$gene_2_left_strand;gene_2_right_strand=$gene_2_right_strand\n";
    #print "$gene_2_left_strand\t$ref:$start..$end($len)\t$gene_2_right_strand\n";
  }
}
