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
    if ($type =~ /mRNA/ or $type =~ /UTR/ or $type =~ /exon/ or $type =~/intron/ or $type=~/intergenic/){
      if ($type !~ /intergenic/){ #intergenic, does not have a parent or a name
      ##this feature's attributes
        my %attr = $loc_feature->attributes;
        my $name_arry_ref = $attr{'parent_id'};
        my $name = ${$name_arry_ref}[0];
        $name =~ s/(\w+)\.\d+?$/$1/;
        ##parent's feature attributes
        my @p_feat = $db->get_feature_by_name($name);
        my $p_feat = shift @p_feat;
        my %p_attr = $p_feat->attributes;
        my $p_notes_arry_ref = $p_attr{'Note'};
        my $p_note = ${$p_notes_arry_ref}[0];
        if ($p_note =~ /transposon/){
          $type = "transposon_$type";
        }
     } 
     if ( ($type eq 'three_prime_UTR' or $type eq 'five_prime_UTR') and (exists $each{$start}{exon} or exists $each{$start}{intron}) ){
        delete $each{$start}{exon} if exists $each{$start}{exon} ;
        delete $each{$start}{intron} if exists $each{$start}{intron} ;
      }elsif ( ($type eq 'transposon_three_prime_UTR' or $type eq 'transposon_five_prime_UTR') and (exists $each{$start}{transposon_exon} or exists $each{$start}{transposon_intron})){
        delete $each{$start}{transposon_exon} if exists $each{$start}{transposon_exon};
        delete $each{$start}{transposon_intron} if exists $each{$start}{transposon_intron};
      }elsif ( ($type eq 'exon' or $type eq 'intron') and (exists $each{$start}{three_prime_UTR} or exists $each{$start}{five_prime_UTR}) ){
	next;
      }elsif ( ($type eq 'transposon_exon' or $type eq 'transposon_intron')  and (exists $each{$start}{transposon_three_prime_UTR} or exists $each{$start}{transposon_five_prime_UTR}) ){
	next;
      }
      $each{$start}{$type}{source}=$source;      
      $each{$start}{$type}{dfs}=$distance_from_start;      
      $each{$start}{$type}{dfe}=$distance_from_end;
      $each{$start}{$type}{insert_type}=$insert_type;
    }
  }
}
#print Dumper \%each;
foreach my $start (keys %each){
  foreach my $insert_feature (keys %{$each{$start}}){
    my $dfs = $each{$start}{$insert_feature}{dfs};
    my $dfe = $each{$start}{$insert_feature}{dfe};
    my $source = $each{$start}{$insert_feature}{source};
    my $insert_type = $each{$start}{$insert_feature}{insert_type};
    push @{$insert_dfs{$insert_feature}} , $dfs; 
    push @{$insert_dfe{$insert_feature}} , $dfe; 
    $features{$insert_feature}++;
    $strains{$source}{insert_feature}{$insert_feature}++;
    $strains{$source}{transposon}{$insert_type}++ if $insert_feature =~ /transposon/;
    $strains{$source}{insert_type}{$insert_type}++ if $insert_feature !~ /mRNA/;
    $inserts{$insert_type}{$insert_feature}++;
    $insert_feature = '';
    $insert_type = '';
  }
}
print "how many mping insertions are there in a specific feature type?\n";
print "feature\tcount\n";
foreach my $feature (sort keys %features){
	my $count = $features{$feature};
        next if $feature =~ /mRNA/;
	print "$feature\t$count\n";
}
print "How many of each type of mping insertion is in each strain?\n";
print "strain\ttype\tcount\n";
foreach my $strain (sort keys %strains){
  foreach my $type (sort keys %{$strains{$strain}{insert_type}}){
	my $count = $strains{$strain}{insert_type}{$type};
	print "$strain\t$type\t$count\n";
  }
  #foreach my $type (sort keys %{$strains{$strain}{transposon}}){
  #	my $count = $strains{$strain}{transposon}{$type};
  #	print "$strain\ttransposon-insertion\t$type\t$count\n";
  #}
  
}
print "how many mping insertions are there in a specific feature type in each strain?\n";
print "strain\tfeature\tcount\n";
foreach my $strain (sort keys %strains){
  foreach my $feature (sort keys %{$strains{$strain}{insert_feature}}){
	my $count = $strains{$strain}{insert_feature}{$feature};
	print "$strain\t$feature\t$count\n";
  }
}
print "is there a preference for what type of insertions is in a specific feature?\n";
print "mping_type\tfeat_type\tcount\n";
foreach my $type (sort keys %inserts){
  foreach my $feat_type  (sort keys %{$inserts{$type}}){
    my $count = $inserts{$type}{$feat_type};
    print "$type\t$feat_type\t$count\n";
  }
}
print "is there a preference for the distance from the start of the feature in which mping is inserted?\n";
print "inserts in features are binned by distance from start of feature\n";
foreach my $type (sort keys %insert_dfs){
	my @dfs = @{$insert_dfs{$type}};
	my @dfe = @{$insert_dfe{$type}};
	#my $partitions = 100;
	my $partitions = [50,100,150,200,250,350,400,450,500,550,1000,2000,5000];
	my $stat_dfs = Statistics::Descriptive::Full->new();
        $stat_dfs->add_data(@dfs);
        my $bins_dfs = $stat_dfs->frequency_distribution_ref($partitions);	

	my $stat_dfe = Statistics::Descriptive::Full->new();
        $stat_dfe->add_data(@dfe);
        my $bins_dfe = $stat_dfe->frequency_distribution_ref($partitions);	

	print "**$type:dfs\n";
	foreach my $bin (sort {$a <=> $b} keys %{$bins_dfs}){
       		print "\t\t$bin\tcount = ${$bins_dfs}{$bin}\n";
   	}
	print "**$type:dfe\n";
	foreach my $bin (sort {$a <=> $b} keys %{$bins_dfe}){
       		print "\t\t$bin\tcount = ${$bins_dfe}{$bin}\n";
   	}
}
