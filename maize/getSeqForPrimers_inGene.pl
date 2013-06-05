#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;

my $sqlite = shift;
my $TEs_file = shift;
die "Please provide a SQLite datafile of a seqfeature db" if !defined $sqlite;

###==> 26_emptySites.plusEvid.sigDifExp.coord.list <==
###1089_ZM_Tourist_289     chr10:86536576..86536893
###11252_ZM_Tourist_256    chr2:171261111..171261395
###11398_ZM_Tourist_293    chr2:177871230..177871461
###11642_ZM_Tourist_26     chr2:189491429..189491607
###12451_ZM_Tourist_205    chr2:215522244..215522529
###13920_ZM_Tourist_299    chr3:30135160..30135464
###14313_ZM_Tourist_289    chr3:57230315..57230632
###16612_ZM_Tourist_4      chr3:211976574..211976938
###17669_ZM_Tourist_245    chr4:17287840..17288127
###20245_ZM_Tourist_405    chr4:198775813..198776127

my %TEs;
open IN, $TEs_file , or die "Can't open $TEs_file\n";
while (my $line = <IN>){
  chomp $line;
  my ($TE, $ref, $start, $end);
  ##don't want any UN or weird chrs, only 1,2,..
  if ($line =~ /\d+:\d+/){
    ($TE, $ref, $start, $end) = $line =~ /(\S+)\s\D*(\d+):(\d+)\.\.(\d+)/; 
    $TEs{$TE}{ref}=$ref;
    $TEs{$TE}{start}=$start;
    $TEs{$TE}{end}=$end;
  }
}
#############################################################################
## make a Bio::DB::SeqFeature::Store object (contains info about  organism in
## SQLite database and all the methods needed to interact with the database
#############################################################################

#Open the sequence database
my $db_obj = Bio::DB::SeqFeature::Store->new(
  -adaptor => 'DBI::SQLite',
  -dsn     => $sqlite
);

print "TE\tTE_coord\tgene_name(strand)(coord)\tseq1(coord)\[\]seq2(coord)\tseq1\[\]seq2\n";
foreach my $TE (keys %TEs){
  my $TE_ref = $TEs{$TE}{ref};
  my $TE_start = $TEs{$TE}{start};
  my $TE_end = $TEs{$TE}{end};
  $TEs{$TE}{gene}{dist} = 100000;## this is bigger than 1kb, genes need to be within 1kb
  my @features_1kb = $db_obj->get_features_by_location($TE_ref,$TE_start - 1000 , $TE_end + 1000);
  my %genes;
  foreach my $feature (@features_1kb) {
    next unless $feature->type =~ /gene/;
    my $gene_id = $feature->primary_id;
    my $gene_name  = $feature->name;
    my $gene_start = $feature->start;
    my $gene_end   = $feature->end;
    my $gene_ref        = $feature->ref;
    my $strand = $feature->strand;
    #my %attr = $feature->attributes;
    #my ($note) = ${$attr{Note}}[0];
    my $gene_info = "$gene_name\t$strand\t$gene_ref:$gene_start..$gene_end";
    my $dist ;
    if ($strand > 0){
      $dist = $gene_start - $TE_end;
    }else {
      $dist = $TE_start - $gene_end;
    }
    if ($dist < $TEs{$TE}{gene}{dist}){
      #$TEs{$TE}{gene}{info}=$gene_info;
      $TEs{$TE}{gene}{name}=$gene_name;
      $TEs{$TE}{gene}{dist}=$dist;
      $TEs{$TE}{gene}{start}=$gene_start;
      $TEs{$TE}{gene}{strand}=$strand;
      $TEs{$TE}{gene}{end}=$gene_end;
      $TEs{$TE}{gene}{ref}=$gene_ref;
      $TEs{$TE}{gene}{id}=$gene_id;
    }## else don't reset it
  }
  my $nearest_gene = $TEs{$TE}{gene}{name} ; 
  my $ng_start = $TEs{$TE}{gene}{start};
  my $ng_end = $TEs{$TE}{gene}{end};
  my $ng_strand =  $TEs{$TE}{gene}{strand};
  my $ng_ref = $TEs{$TE}{gene}{ref} ; 
  my $ng_id = $TEs{$TE}{gene}{id} ; 
  next if !defined $ng_id;
  my @f = $db_obj->get_features_by_location($ng_ref,$ng_start,$ng_end);
  foreach my $f (@f){
    next if $f->type !~ /gene/;
    my $f_id  = $f->primary_id;
    next if $ng_id ne $f_id ;
    my $f_name  = $f->name;
    my $f_ref  = $f->ref;
    my $f_start = $f->start;
    my $f_end   = $f->end;
    my $f_strand = $f->strand;
    my ($seq1_type,$seq2_type,$seq1_start,$seq1_end,$seq2_start,$seq2_end);
    my $upstream_of_TE_bp = 100;
    my $start_of_gene_bp = 300;
    ##start will always be less than end
    if ($f_strand > 0){
      ## seq1 => TE
      ## seq2 => gene
      ## ===TE==== >>>>>gene>>>>>
      $seq1_type = "upstream_of_TE";
      $seq2_type = "start_of_gene";
      $seq2_start = $f_start;
      $seq2_end   = $f_start + $start_of_gene_bp;
      $seq1_start   = $TE_start - $upstream_of_TE_bp;
      $seq1_end     = $TE_start - 1;
    }else {
      ## seq1 => gene 
      ## seq2 => TE
      ## <<<<<gene<<<<<  ====TE====
      $seq1_type = "start_of_gene";
      $seq2_type = "upstream_of_TE";
      $seq1_start = $f_end - $start_of_gene_bp ; 
      $seq1_end   = $f_end;
      $seq2_start   = $TE_end + 1;
      $seq2_end     = $TE_end + $upstream_of_TE_bp;
    }
    my $strand = $f_strand > 0 ? '+' : '-' ; 
    my $seq1 =$db_obj->fetch_sequence(-seq_id=>$TE_ref,-start=>$seq1_start,-end=>$seq1_end);
    my $seq2 =$db_obj->fetch_sequence(-seq_id=>$TE_ref,-start=>$seq2_start,-end=>$seq2_end);

    print "$TE\t$TE_ref:$TE_start..$TE_end\t$f_name($strand)($f_ref:$f_start..$f_end)\t$seq1_type($f_ref:$seq1_start..$seq1_end)\[\]$seq2_type($f_ref:$seq2_start..$seq2_end)\t$seq1\[\]$seq2\n";
 
  }
}
