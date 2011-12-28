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
my %features;
my %strains;
my %transposons;
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
    my $dfs = '';
    foreach my $loc_feature (@loc_features){
    	my $type = $loc_feature->type;
	my $f_start = $loc_feature->start;
	my $f_end = $loc_feature->end;
	my $f_strand = $loc_feature->strand;
	my $f_len = $loc_feature->length;
        if ($f_strand eq '+'){
		$distance_from_start = ($start - $f_start) 
	}else{
	  	$distance_from_start = ($f_end - $start); 
        }
##the 9th column Note=
        my $transponCheck = 1;
        if ($transponCheck){
          my %attr = $loc_feature->attributes;
            if ($type =~ /gene/){
               my $notes_arry_ref = $attr{'Note'};
               my $note = ${$notes_arry_ref}[0];
               my $name = ${$attr{'load_id'}}[0];
               $transposons{$name}=0 if $note =~ /transposon/;
	       $insert_feature = 'transposon';
               $dfs = $distance_from_start;
    	       #push @{$insert_dfs{'transposon'}} , $distance_from_start; 
               last if $note =~ /transposon/;
             }elsif($type =~ /mRNA/ or $type =~ /UTR/ or $type =~ /exon/){
               my $name_arry_ref = $attr{'parent_id'};
               my $name = ${$name_arry_ref}[0];
               $name =~ s/(\w+)\.\d+?$/$1/;
	       $insert_feature = 'transposon';
               $dfs = $distance_from_start;
	       
    	       #push @{$insert_dfs{'transposon'}} , $distance_from_start; 
               last if exists $transposons{$name};
             }
        }
## end transposon check
        if ($type =~ /mRNA/){
    	  push @{$insert_dfs{'mRNA'}} , $distance_from_start; 
	}
	if ($type =~ /three_prime_UTR/){
    	  push @{$insert_dfs{'three_prime_UTR'}} , $distance_from_start; 
        }
	if ($type =~ /five_prime_UTR/){
    	  push @{$insert_dfs{'five_prime_UTR'}} , $distance_from_start; 
        }
	if ($type =~ /exon/){
	  $insert_feature = 'exon';
          $dfs = $distance_from_start;
        }elsif ($type eq 'intron'){
	  $insert_feature = $type;
          $dfs = $distance_from_start;
        }elsif($type eq 'intergenic'){
	  $insert_feature = $type;
          $dfs = $distance_from_start;
        }
    } 
    push @{$insert_dfs{$insert_feature}} , $dfs; 
    $inserts{$insert_type}{$insert_feature}++;
    $features{$insert_feature}++;
    $strains{$source}{insert_feature}{$insert_feature}++;
    $strains{$source}{insert_type}{$insert_type}++;
    $strains{$source}{transposon}{$insert_type}++ if $insert_feature eq 'transposon';
    #print "$source $ref:$start..$end $insert_type $avg_flankers $spanners $insert_feature\n"; 

}
print "how many mping insertions are there in a specific feature type?\n";
print "feature\tcount\n";
foreach my $feature (sort keys %features){
	my $count = $features{$feature};
	print "$feature\t$count\n";
}
print "How many of each type of mping insertion is in each strain?\n";
print "strain\ttype\tcount\n";
foreach my $strain (sort keys %strains){
  foreach my $type (sort keys %{$strains{$strain}{insert_type}}){
	my $count = $strains{$strain}{insert_type}{$type};
	print "$strain\t$type\t$count\n";
  }
  foreach my $type (sort keys %{$strains{$strain}{transposon}}){
	my $count = $strains{$strain}{transposon}{$type};
	print "$strain\ttransposon-insertion\t$type\t$count\n";
  }
  
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
print "inserts in features are binned by relative distance from start of feature ((f.start-insert.pos)/f.len)*100\n";
foreach my $type (sort keys %insert_dfs){
	my @dfs = @{$insert_dfs{$type}};
	#my $partitions = 100;
	my $partitions = [50,100,150,200,250,350,400,450,500,550,1000,2000,5000];
	my $stat_dfs = Statistics::Descriptive::Full->new();
        $stat_dfs->add_data(@dfs);
        my $bins_dfs = $stat_dfs->frequency_distribution_ref($partitions);	

	print "**$type\n";
	foreach my $bin (sort {$a <=> $b} keys %{$bins_dfs}){
       		print "\t\t$bin\tcount = ${$bins_dfs}{$bin}\n";
   	}
}
