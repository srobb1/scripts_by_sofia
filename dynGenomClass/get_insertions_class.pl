#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;
use Statistics::Descriptive;  
use File::Spec;

my $sqlite = shift;
die "Please provide a SQLite datafile of a seqfeature db" if !defined $sqlite;
$sqlite = File::Spec->rel2abs($sqlite);

# Open the sequence database
  my $db = Bio::DB::SeqFeature::Store->new(-adaptor => 'DBI::SQLite',
                                           -dsn     => $sqlite);

# queries DB for all features of type 'transposable_element_insertion_site'
# transposable_element_attribute'
#my $query_type = 'transposable_element_insertion_site';
my $query_type = 'transposable_element_attribute';
my %inserts;
my %insert_dfs;
my %insert_dfe;
my %features;
my %strains;
my %each;
my @features_type = $db->get_features_by_type($query_type);
foreach my $feature (@features_type) {
  my $source = $feature->source;	
  my $name  = $feature->name;
  my $start = $feature->start;
  my $end   = $feature->end;
  my $ref        = $feature->ref;
  #ID=Chr8.9323600.spanners;avg_flankers=10.5;spanners=12;type='heterozygous';
  my %attr_hash = $feature->attributes;
  my $insert_type =  ${$attr_hash{'type'}}[0];
  $insert_type = 'homozygous' if $insert_type =~ /homo/; 
  $insert_type = 'heterozygous' if $insert_type =~ /het/; 
  $insert_type = 'somatic' if $insert_type =~ /new/; 
  my $avg_flankers =  ${$attr_hash{'avg_flankers'}}[0];
  my $spanners =  ${$attr_hash{'spanners'}}[0];

  #what feature is the insertion inside of
  my @loc_features = $db->get_features_by_location(-seq_id=>$ref,-start=>$start,-end=>$end);
  my $insert_feature = '';
  my $distance_from_start = '';
  my $distance_from_end = '';
  my $dfs = '';
  my $dfe = '';
  foreach my $loc_feature (@loc_features){
    my $type = $loc_feature->type;
    $type =~ s/:.+$//;
    my $f_start = $loc_feature->start;
    my $f_end = $loc_feature->end;
    my $f_strand = $loc_feature->strand;
    my $f_len = $loc_feature->length;
    if ($f_strand eq '+'){
      $distance_from_start = ($start - $f_start); 
      $distance_from_end = ($f_end - $start); 
    }else{
      $distance_from_start = ($f_end - $start); 
      $distance_from_end = ($start - $f_start); 
    }
    my ($gene2right, $gene2left,$f_name,$f_note) = qw (n/a n/a n/a n/a);
    if ($type =~ /mRNA/ or $type =~ /UTR/ or $type =~ /exon/ or $type =~/intron/ or $type=~/intergenic/){
      if ($type !~ /intergenic/){ #intergenic, does not have a parent or a name
      ##this feature's attributes
        my %attr = $loc_feature->attributes;
        my $name_arry_ref = $attr{'parent_id'};
        $f_name = ${$name_arry_ref}[0];
        $f_name =~ s/(\w+)\.\d+?$/$1/;
        ##parent's feature attributes
        my @p_feat = $db->get_feature_by_name($f_name);
        my $p_feat = shift @p_feat;
        my %p_attr = $p_feat->attributes;
        my $p_notes_arry_ref = $p_attr{'Note'};
        my $p_note = ${$p_notes_arry_ref}[0];
        $f_note=$p_note;
        if ($p_note =~ /transposon/){
          $type = "transposon_$type";
        }
     }else {
       ##if intergenic, what is strand of surrounding genes
       ##upstream_gene,downstream_gene  gene_2_right_strand, gene_2_left_strand
       my %attr = $loc_feature->attributes;
       my $gene2right_ref = $attr{'upstream_gene'};
       my $gene2left_ref = $attr{'downstream_gene'};
       my $upstream = ${$gene2right_ref}[0];
       my $downstream = ${$gene2left_ref}[0];
       $upstream =~ s/.+\((\+|\-)\)$/$1/;
       $downstream =~ s/.+\((\+|\-)\)$/$1/;
       $gene2right = $upstream eq '+' ? 1 : -1;
       $gene2left = $downstream eq '+' ? 1 : -1;
     } 
     if ( ($type eq 'three_prime_UTR' or $type eq 'five_prime_UTR') and (exists $each{$ref}{$start}{exon} or exists $each{$ref}{$start}{intron}) ){
        delete $each{$ref}{$start}{exon} if exists $each{$ref}{$start}{exon} ;
        delete $each{$ref}{$start}{intron} if exists $each{$ref}{$start}{intron} ;
      }elsif ( ($type eq 'transposon_three_prime_UTR' or $type eq 'transposon_five_prime_UTR') and (exists $each{$ref}{$start}{transposon_exon} or exists $each{$ref}{$start}{transposon_intron})){
        delete $each{$ref}{$start}{transposon_exon} if exists $each{$ref}{$start}{transposon_exon};
        delete $each{$ref}{$start}{transposon_intron} if exists $each{$ref}{$start}{transposon_intron};
      }elsif ( ($type eq 'exon' or $type eq 'intron') and (exists $each{$ref}{$start}{three_prime_UTR} or exists $each{$ref}{$start}{five_prime_UTR}) ){
	next;
      }elsif ( ($type eq 'transposon_exon' or $type eq 'transposon_intron')  and (exists $each{$ref}{$start}{transposon_three_prime_UTR} or exists $each{$ref}{$start}{transposon_five_prime_UTR}) ){
	next;
      }
      $each{$ref}{$start}{$type}{spanners}=$spanners;      
      $each{$ref}{$start}{$type}{flankers}=$avg_flankers;      
      $each{$ref}{$start}{$type}{f_end}=$f_end;      
      $each{$ref}{$start}{$type}{f_start}=$f_start;      
      $each{$ref}{$start}{$type}{f_strand}=$f_strand;      
      $each{$ref}{$start}{$type}{f_name}=$f_name;      
      $each{$ref}{$start}{$type}{f_note}=$f_note;      
      $each{$ref}{$start}{$type}{source}=$source;      
      $each{$ref}{$start}{$type}{dfs}=$distance_from_start;      
      $each{$ref}{$start}{$type}{dfe}=$distance_from_end;
      $each{$ref}{$start}{$type}{feat_len}=$f_len;
      $each{$ref}{$start}{$type}{insert_type}=$insert_type;
      if ($type eq 'intergenic'){
        $each{$ref}{$start}{$type}{gene2right}=$gene2right;
        $each{$ref}{$start}{$type}{gene2left}=$gene2left;
   
      }
    }
  }
}


print "source\tref\tstart\tinsert_type\tspanners\tavg.flankers\tinsert_feature\tf_start\tf_end\tf_strand\tf_name\tf_note\n";
foreach my $ref (keys %each){
  foreach my $start (keys %{$each{$ref}}){
    foreach my $insert_feature (keys %{$each{$ref}{$start}}){
      my $spanners = $each{$ref}{$start}{$insert_feature}{spanners};      
      my $flankers = $each{$ref}{$start}{$insert_feature}{flankers};      
      my $dfs = $each{$ref}{$start}{$insert_feature}{dfs};
      my $dfe = $each{$ref}{$start}{$insert_feature}{dfe};
      my $f_start = $each{$ref}{$start}{$insert_feature}{f_start};
      my $f_end = $each{$ref}{$start}{$insert_feature}{f_end};
      my $f_strand= $each{$ref}{$start}{$insert_feature}{f_strand};
      my $source = $each{$ref}{$start}{$insert_feature}{source};
      my $f_name = $each{$ref}{$start}{$insert_feature}{f_name};
      my $f_note = $each{$ref}{$start}{$insert_feature}{f_note};
      my $insert_type = $each{$ref}{$start}{$insert_feature}{insert_type};
      my $intergenic_g2r = $each{$ref}{$start}{$insert_feature}{gene2right};
      my $intergenic_g2l = $each{$ref}{$start}{$insert_feature}{gene2left};
      if ($insert_feature eq 'intergenic'){
	$f_strand = "$intergenic_g2l/$intergenic_g2r";
      }
      print "$source\t$ref\t$start\t$insert_type\t$spanners\t$flankers\t$insert_feature\t$f_start\t$f_end\t$f_strand\t$f_name\t$f_note\n";
    }
  }
}
