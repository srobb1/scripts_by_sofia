#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);

# Open the sequence database
my $db_file = shift;
my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::SQLite',
                                           #-dsn     => '/home_stajichlab/robb/rice/database/MSUr7.mping.insertions.sqlite');
                                           -dsn     => $db_file);
my %ref;
my $ref_total = 0;
my %transposons;
my %feat;
#my %strand;

my @features_type = $db->get_features_by_type('chromosome');
foreach my $feature (@features_type) {
  my $ref = $feature->ref;
  next if $ref =~/Un/ or $ref =~ /Sy/;
  my $len = $feature->length;
  $ref{$ref}= $len ;
  $ref_total+=$len;
}
foreach my $ref (keys %ref){
    for (my $i = 0; $i < $ref{$ref} + 1 ; $i++){
      ## set every bp to 1
      ${$feat{intergenic}{$ref}}[$i] = 1;
      ## to keep track of gene strand
      #${$strand{$ref}}[$i] = 0;
    }
    ${$feat{intergenic}{$ref}}[0] = 0;
}
print "ref: $ref_total\n";
foreach my $feature ($db->get_features_by_type('gene')) {
  my $ref = $feature->ref; 
  next if !exists $ref{$ref}; 
# if (!exists $feat{intergenic}{$ref}){
 #   for (my $i = 0; $i < $ref{$ref} + 1 ; $i++){
 #     ## set every bp to 1
 #     ${$feat{intergenic}{$ref}}[$i] = 1;
 #     ## to keep track of gene strand
 #     ${$strand{$ref}}[$i] = 0;
 #     ${$feat{intergenic}{$ref}}[0] = 0;
 #   }
 # }
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
    #${$strand{$ref}}[$i] = $f_strand; 
    ## where a gene exists change intergenic index to 0
    my $strand = '*';
    if ($f_strand < 0){
      $strand = '-';
    }elsif ($f_strand > 0){
      $strand = '+';
    }
    ${$feat{intergenic}{$ref}}[$i] = $strand; 
  }
}  

  my $feat_total = 0;
  foreach my $ref (sort keys %{$feat{intergenic}}){
    pop @{$feat{intergenic}{$ref}}; #remove bp index 0, doesnt really exisit
    my $str = join '' , @{$feat{intergenic}{$ref}};
    $str =~ s/\D//g;
    my $feat_sum = length $str;
    #my $feat_sum = sum @{$feat{intergenic}{$ref}}; 
    $feat_total += $feat_sum;
  }
  print "intergenic: $feat_total\n";
my %intergenic;
foreach my $ref (sort keys %ref){
  #my $last_value = 0;
  #my $start;
  #my $end;
  #my $gene_2_right_strand;
  #my $gene_2_left_strand;

#  for (my $i=0; $i < $ref{$ref} + 1 ; $i++){
#    my $value = ${$feat{intergenic}{$ref}}[$i];
#      if ($last_value == 0 and $value == 1){
#        ##first position of intergenic
#        $start = $i;
#        if ($i == 0){
#          $gene_2_left_strand = 0;
#        }else {
#          $gene_2_left_strand = $strand{$ref}[$i-1]
#        }
#      }elsif($last_value == 1 and $value == 0) {
#        ##last position of intergenic
#        $end = $i-1;
#        $gene_2_right_strand = $strand{$ref}[$i];
#        $intergenic{$ref}{$start}{end} = $end;
#        $intergenic{$ref}{$start}{left_strand} = $gene_2_left_strand;
#        $intergenic{$ref}{$start}{right_strand} = $gene_2_right_strand;
#        ## reset variables
#        $gene_2_right_strand = '';
#        $gene_2_left_strand = '';
#	$start = '';
#        $end = '';
#      }
#      $last_value = $value;
#  }

  my $intergenic = join '' , @{$feat{intergenic}{$ref}};
  while ($intergenic =~ /(\-|\+)*(1+)(\-|\+)*/g){
    ##postion to the right of match
    ##0-base notation
    ##need to calcuate for 1-base notation
    my $end = pos $intergenic;
    my $len = length $2;
    my $start = $end - $len + 1;
    my $g2l = $1;
    my $g2r = $3;
    my $g2l_strand = 0;
    my $g2r_strand = 0;
    if (defined $g2l){
      my @g2l = split '',$g2l;
      $g2l_strand = pop @g2l;
      if ($g2l_strand eq '-'){
        $g2l_strand = -1;
      }elsif ($g2l_strand eq '+'){
        $g2l_strand = 1;
      }
    }
    if (defined $g2r){
      my @g2r = split '',$g2r;
      $g2r_strand = shift @g2r;
      if ($g2r_strand eq '-'){
        $g2r_strand = -1;
      }elsif ($g2r_strand eq '+'){
        $g2r_strand = 1;
      }
    }
    $intergenic{$ref}{$start}{end} = $end;
    $intergenic{$ref}{$start}{left_strand} = $g2l_strand;
    $intergenic{$ref}{$start}{right_strand} = $g2r_strand;
  }
}
open GFF, ">intergenic.gff" or die "Can't open intergenic.gff for writing\n";
my $promoters;
my $threeprime;
foreach my $ref (sort keys %intergenic){
  foreach my $start (sort {$a <=> $b} keys %{$intergenic{$ref}}){
    my $end = $intergenic{$ref}{$start}{end};
    my $gene_2_left_strand = $intergenic{$ref}{$start}{left_strand};
    my $gene_2_right_strand = $intergenic{$ref}{$start}{right_strand};
    my $len = $end - $start + 1;
    print GFF "$ref\t.\tintergenic\t$start\t$end\t.\t.\t.\tID=$ref.intergenic.$start.$end;gene_2_left_strand=$gene_2_left_strand;gene_2_right_strand=$gene_2_right_strand\n";
    if ($gene_2_right_strand > 0 ){
      $promoters++;
    }if ($gene_2_left_strand < 0){
      $promoters++;
    }
    if ($gene_2_right_strand < 0 ){
      $threeprime++;
    }if ($gene_2_left_strand > 0){
      $threeprime++;
    }
    #print "$gene_2_left_strand\t$ref:$start..$end($len)\t$gene_2_right_strand\n";
  }
}
print "promoter count = $promoters, promoter bp = ",$promoters * 500 ,"\n";
print "threeprime count = $threeprime, threeprime bp = ",$threeprime * 500 ,"\n";
