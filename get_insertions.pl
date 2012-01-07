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
    my ($gene2right, $gene2left);
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
      $each{$ref}{$start}{$type}{source}=$source;      
      $each{$ref}{$start}{$type}{dfs}=$distance_from_start;      
      $each{$ref}{$start}{$type}{dfe}=$distance_from_end;
      $each{$ref}{$start}{$type}{feat_len}=$f_len;
      $each{$ref}{$start}{$type}{insert_type}=$insert_type;
      $each{$ref}{$start}{$type}{f_strand}=$f_strand;
      if ($type eq 'intergenic'){
        $each{$ref}{$start}{$type}{gene2right}=$gene2right;
        $each{$ref}{$start}{$type}{gene2left}=$gene2left;
   
      }
    }
  }
}
#print Dumper \%each;
foreach my $ref (keys %each){
  foreach my $start (keys %{$each{$ref}}){
    foreach my $insert_feature (keys %{$each{$ref}{$start}}){
      my $dfs = $each{$ref}{$start}{$insert_feature}{dfs};
      my $dfe = $each{$ref}{$start}{$insert_feature}{dfe};
      my $source = $each{$ref}{$start}{$insert_feature}{source};
      my $insert_type = $each{$ref}{$start}{$insert_feature}{insert_type};
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
my $toPrint = "souce\tref\tinsert_type\tinsert_feature\tf_len\tdfs\tdfe\trel_pos\tf_strand\n";	
my %intergenic;
foreach my $ref (sort keys %each){
  foreach my $start (sort keys %{$each{$ref}}){
    foreach my $insert_feature (sort keys %{$each{$ref}{$start}}){
      my $dfs = $each{$ref}{$start}{$insert_feature}{dfs};
      my $dfe = $each{$ref}{$start}{$insert_feature}{dfe};
      my $source = $each{$ref}{$start}{$insert_feature}{source};
      my $insert_type = $each{$ref}{$start}{$insert_feature}{insert_type};
      my $f_len = $each{$ref}{$start}{$insert_feature}{feat_len};
      my $f_strand = $each{$ref}{$start}{$insert_feature}{f_strand};
      my $rel_pos = $dfs / $f_len ;
      my ($g2r, $g2l) = '';
      my $orientation;
      if ($insert_feature eq 'intergenic'){
        $g2r = $each{$ref}{$start}{$insert_feature}{gene2right};
        $g2l = $each{$ref}{$start}{$insert_feature}{gene2left};
        if ($g2l == $g2r){
          $orientation = 1;
        }elsif ($g2l > 0 and $g2r < 0){
          $orientation = 0;
        }elsif($g2l < 0 and $g2r > 0) {
          $orientation = 2;
        }else {
	  $orientation = "?**$g2l**$g2r**?";
	}
	$f_strand = $orientation;
      }
      $toPrint .= "$source\t$ref\t$insert_type\t$insert_feature\t$f_len\t$dfs\t$dfe\t$rel_pos\t$f_strand\n";	
      ##need to characterized the intergenic insertion type:
      ## ---->(+) V (+)----->  tandem: 1 only inserted into 1 5' gene region
      ## <----(-) V (-)<-----  tandem: 1
      ## ---->(+) V (-)<-----  convergent: 0
      ## <----(-) V (+)----->  diverent : 2
      ##also i need to look at insertions that are not in the 5kb regions before or after a gene     

      ##look at distance from 5' end in tandem intergenic regions, therefor only 1 5'end
      ##also want to plot the distances from the gene boundries if the intergenic region is inter<5000, 11000<inter>5000, inter>11000


      if ($insert_feature eq 'intergenic' and $f_len > 5000){
	## distance from upstream gene
        push @{$intergenic{all}{gt_5000}{dug}} , $dfe if $g2r > 0;#g2r == +1
        push @{$intergenic{all}{gt_5000}{dug}} , $dfs if $g2l < 0;#g2l == -1
        push @{$intergenic{all}{gt_5000}{ddg}} , $dfe if $g2r < 0;#g2r == -1
        push @{$intergenic{all}{gt_5000}{ddg}} , $dfs if $g2l > 0;#g2l == +1
        push @{$intergenic{all}{gt_5000}{dfe}} , $dfe;
        push @{$intergenic{all}{gt_5000}{dfs}} , $dfs;
        push @{$intergenic{$insert_type}{gt_5000}{dug}} , $dfe if $g2r > 0;#g2r == +1
        push @{$intergenic{$insert_type}{gt_5000}{dug}} , $dfs if $g2l < 0;#g2l == -1
        push @{$intergenic{$insert_type}{gt_5000}{ddg}} , $dfe if $g2r < 0;#g2r == -1
        push @{$intergenic{$insert_type}{gt_5000}{ddg}} , $dfs if $g2l > 0;#g2l == +1
        push @{$intergenic{$insert_type}{gt_5000}{dfe}} , $dfe;
        push @{$intergenic{$insert_type}{gt_5000}{dfs}} , $dfs;

        if ($orientation == 1){
          push @{$intergenic{all}{gt_5000}{df5_tandem}}, $dfs if $g2l > 0;
          push @{$intergenic{all}{gt_5000}{df5_tandem}}, $dfe if $g2r > 0;
          push @{$intergenic{all}{gt_5000}{df3_tandem}}, $dfs if $g2l < 0;
          push @{$intergenic{all}{gt_5000}{df3_tandem}}, $dfe if $g2r > 0;
          push @{$intergenic{$insert_type}{gt_5000}{df5_tandem}}, $dfs if $g2l > 0;
          push @{$intergenic{$insert_type}{gt_5000}{df5_tandem}}, $dfe if $g2r > 0;
          push @{$intergenic{$insert_type}{gt_5000}{df3_tandem}}, $dfs if $g2l < 0;
          push @{$intergenic{$insert_type}{gt_5000}{df3_tandem}}, $dfe if $g2r > 0;
        }elsif($orientation == 2){
          push @{$intergenic{all}{gt_5000}{dfs_divergent}}, $dfs;
          push @{$intergenic{all}{gt_5000}{dfe_divergent}}, $dfe;
          push @{$intergenic{$insert_type}{gt_5000}{dfs_divergent}}, $dfs;
          push @{$intergenic{$insert_type}{gt_5000}{dfe_divergent}}, $dfe;
        }elsif ($orientation == 0){
          push @{$intergenic{all}{gt_5000}{dfs_convergent}}, $dfs;
          push @{$intergenic{all}{gt_5000}{dfe_convergent}}, $dfe;
          push @{$intergenic{$insert_type}{gt_5000}{dfs_convergent}}, $dfs;
          push @{$intergenic{$insert_type}{gt_5000}{dfe_convergent}}, $dfe;
        }

      }elsif ($insert_feature eq 'intergenic' and $f_len < 5000){
        push @{$intergenic{all}{lt_5000}{dug}} , $dfe if $g2r > 0;#g2r == +1
        push @{$intergenic{all}{lt_5000}{dug}} , $dfs if $g2l < 0;#g2l == -1
        push @{$intergenic{all}{lt_5000}{ddg}} , $dfe if $g2r < 0;#g2r == -1
        push @{$intergenic{all}{lt_5000}{ddg}} , $dfs if $g2l > 0;#g2l == +1
        push @{$intergenic{all}{lt_5000}{dfe}} , $dfe;
        push @{$intergenic{all}{lt_5000}{dfs}} , $dfs;
        push @{$intergenic{$insert_type}{lt_5000}{dug}} , $dfe if $g2r > 0;#g2r == +1
        push @{$intergenic{$insert_type}{lt_5000}{dug}} , $dfs if $g2l < 0;#g2l == -1
        push @{$intergenic{$insert_type}{lt_5000}{ddg}} , $dfe if $g2r < 0;#g2r == -1
        push @{$intergenic{$insert_type}{lt_5000}{ddg}} , $dfs if $g2l > 0;#g2l == +1
        push @{$intergenic{$insert_type}{lt_5000}{dfe}} , $dfe;
        push @{$intergenic{$insert_type}{lt_5000}{dfs}} , $dfs;
        if ($orientation == 1){
          push @{$intergenic{all}{lt_5000}{df5_tandem}}, $dfs if $g2l > 0;
          push @{$intergenic{all}{lt_5000}{df5_tandem}}, $dfe if $g2r > 0;
          push @{$intergenic{all}{lt_5000}{df3_tandem}}, $dfs if $g2l < 0;
          push @{$intergenic{all}{lt_5000}{df3_tandem}}, $dfe if $g2r > 0;
          push @{$intergenic{$insert_type}{lt_5000}{df5_tandem}}, $dfs if $g2l > 0;
          push @{$intergenic{$insert_type}{lt_5000}{df5_tandem}}, $dfe if $g2r > 0;
          push @{$intergenic{$insert_type}{lt_5000}{df3_tandem}}, $dfs if $g2l < 0;
          push @{$intergenic{$insert_type}{lt_5000}{df3_tandem}}, $dfe if $g2r > 0;
        }elsif($orientation == 2){
          push @{$intergenic{all}{lt_5000}{dfs_divergent}}, $dfs;
          push @{$intergenic{all}{lt_5000}{dfe_divergent}}, $dfe;
          push @{$intergenic{$insert_type}{lt_5000}{dfs_divergent}}, $dfs;
          push @{$intergenic{$insert_type}{lt_5000}{dfe_divergent}}, $dfe;
        }elsif ($orientation == 0){
          push @{$intergenic{all}{lt_5000}{dfs_convergent}}, $dfs;
          push @{$intergenic{all}{lt_5000}{dfe_convergent}}, $dfe;
          push @{$intergenic{$insert_type}{lt_5000}{dfs_convergent}}, $dfs;
          push @{$intergenic{$insert_type}{lt_5000}{dfe_convergent}}, $dfe;
        }
      }
    }
  }
}
##distance from up and downstream genes
my %distance;
foreach my $i_type ( sort keys %intergenic ) {
  my @dug        = @{ $intergenic{$i_type}{gt_5000}{dug} };
  my $partitions = [ 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 5500, 6000, 6500, 50000, 100000 ];
  my $desc =  "intergenic gt 5000; distance from 5' end of gene:$i_type";
  my $bins = getBins ($intergenic{$i_type}{gt_5000}{dug},$partitions);   
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic gt 5000; distance from 3' end of gene:$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{ddg},$partitions);   
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from 5' end of gene:$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{dug},$partitions);   
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from 3' end of gene:$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{ddg},$partitions);   
  add2DistanceHash ($bins,$desc);
  

  $desc =  "intergenic gt 5000; distance from 5' end of gene (tandem --> -->):$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{ddg},$partitions);   
  $bins = getBins ($intergenic{$i_type}{gt_5000}{df5_tandem},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic gt 5000; distance from 3' end of gene (tandem --> -->):$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{df3_tandem},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from 5' end of gene (tandem --> -->):$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{df5_tandem},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from 3' end of gene (tandem --> -->):$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{df3_tandem},$partitions);
  add2DistanceHash ($bins,$desc);


  $desc =  "intergenic gt 5000; distance from start of intergenic region:$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{dfs},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic gt 5000; distance from end of intergenic region:$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{dfe},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from start of intergenic region:$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{dfs},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from end of intergenic region:$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{dfe},$partitions);
  add2DistanceHash ($bins,$desc);

  $desc =  "intergenic gt 5000; distance from start of intergenic region (convergent --> <--):$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{dfs_convergent},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic gt 5000; distance from end of intergenic region (convergent --> <--):$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{dfe_convergent},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from start of intergenic region (convergent --> <--):$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{dfs_convergent},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic lt 5000; distance from end of intergenic region (convergent --> <--):$i_type";
  $bins = getBins ($intergenic{$i_type}{lt_5000}{dfe_convergent},$partitions);
  add2DistanceHash ($bins,$desc);
  
  $desc =  "intergenic gt 5000; distance from start of intergenic region (divergent <-- -->):$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{dfs_divergent},$partitions);
  add2DistanceHash ($bins,$desc);
  $desc =  "intergenic gt 5000; distance from end of intergenic region (divergent <-- -->):$i_type";
  $bins = getBins ($intergenic{$i_type}{gt_5000}{dfe_divergent},$partitions);
  add2DistanceHash ($bins,$desc);
  if (defined $intergenic{$i_type}{lt_5000}{dfs_divergent}){
    $desc =  "intergenic lt 5000; distance from start of intergenic region (divergent <-- -->):$i_type";
    $bins = getBins ($intergenic{$i_type}{lt_5000}{dfs_divergent},$partitions);
    add2DistanceHash ($bins,$desc);
  }
  if (defined $intergenic{$i_type}{lt_5000}{dfe_divergent}){
    $desc =  "intergenic lt 5000; distance from end of intergenic region (divergent <-- -->):$i_type";
    $bins = getBins ($intergenic{$i_type}{lt_5000}{dfe_divergent},$partitions);
    add2DistanceHash ($bins,$desc);
  }
}
my @array;
push @{$array[0]} ,  sort keys %distance;
foreach my $desc (sort keys %distance){
    my $i = 1;
    foreach my $bin (sort {$a <=> $b} keys %{$distance{$desc}}){
      $array[$i][0] = $bin;
      my $count = $distance{$desc}{$bin};
      push @{$array[$i]},$count;
      $i++;
    }
}

print "\t";
for (my $i=0 ; $i<@array ; $i++){
  print join ("\t", @{$array[$i]}), "\n"; 
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


print $toPrint;


sub getBins {
  my $array_ref = shift;
  my $partitions_ref = shift;
  my @array        = @{$array_ref}; 
  my $partitions = $partitions_ref;
  my $stat   = Statistics::Descriptive::Full->new();
  $stat->add_data(@array);
  my $bins = $stat->frequency_distribution_ref($partitions);
  return $bins;
}

sub add2DistanceHash { 
  my $bins = shift;
  my $desc = shift;
  foreach my $bin ( sort { $a <=> $b } keys %{$bins} ) {
    my $count = ${$bins}{$bin};
    $distance{$desc}{$bin}=$count;
  }
}
