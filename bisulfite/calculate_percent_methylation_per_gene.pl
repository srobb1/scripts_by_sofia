#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;

my $sqlite = shift;
my $dir    = "/rhome/robb/project/bisulfide/output_all/split_methyl";    #shift

die "Please provide a SQLite datafile of a seqfeature db" if !defined $sqlite;

#############################################################################
## make a Bio::DB::SeqFeature::Store object (contains info about  organism in
## SQLite database and all the methods needed to interact with the database
#############################################################################

#Open the sequence database
my $db_obj =
  Bio::DB::SeqFeature::Store->new( -adaptor => 'DBI::SQLite',
                                   -dsn     => $sqlite );

print "gene\tref:start..end\tstrand\tstrain\t%CpG\t%CHH\t%CHG\tnote\n";
## get_feature_by_type can get any type from your original GFF, anything in Col3
my @features = $db_obj->get_features_by_type('gene');
foreach my $f (@features) {
  my $ref    = $f->ref;
  my $start  = $f->start;
  my $end    = $f->end;
  my $strand = $f->strand;
  my %attr   = $f->attributes;
  my $note   = ${ $attr{'Note'} }[0];
  my $name   = $f->name;

  foreach my $strain (qw (A119 NB RIL16 RIL25 RIL43 RIL60)) {
    my %results;
    foreach my $type (qw (CHH CpG CHG)) {
      my %methy;
      open IN, "$dir/${strain}_${ref}_${type}.split.txt" or die "Can't open $dir/${strain}_${ref}_${type}.split.txt"; #NB_Chr9_CpG.split.txt
      while ( my $line = <IN> ) {
        chomp $line;
        my ( $pos, $methyl_yes, $methyl_no ) = split /\t/, $line;
        next unless $pos >= $start;
        next unless $pos <= $end;
        last if $pos > $end;
        $methyl_yes = defined $methyl_yes ? $methyl_yes : 0;
        $methyl_no  = defined $methyl_no  ? $methyl_no  : 0;
        $methy{$pos}{yes} += $methyl_yes;
        $methy{$pos}{no}  += $methyl_no;
      }
      my $region_yes = 0;
      my $region_no  = 0;
      foreach my $pos ( sort { $a <=> $b } keys %methy ) {
        my ( $methy_yes, $methy_no ) = ( 0, 0 );
        if ( exists $methy{$pos}{yes} ) {
          $methy_yes = $methy{$pos}{yes};
          $region_yes += $methy_yes;
        }
        if ( exists $methy{$pos}{no} ) {
          $methy_no = $methy{$pos}{no};
          $region_no += $methy_no;
        }

        #my $total = $methy_yes + $methy_no;
        #my $percent = ( $methy_yes / ($total) ) * 100;

        #my $pretty_percent = sprintf( '%.1f', $percent );
      }
      my $region_total          = $region_yes + $region_no;
      my $region_percent;
      my $pretty_region_percent; 
      if ($region_total == 0 or !defined $region_total){
        #next;
        #print "$name $strain $ref:$start..$end  rtotal:$region_total r_yes:$region_yes r_no:$region_no\n";
        $pretty_region_percent = 'n/a'; 
      }else{
        $region_percent        = ( $region_yes / $region_total ) * 100;
        $pretty_region_percent = sprintf( '%.1f', $region_percent );
      }

      #$results{$type}{region_yes}=$region_yes;
      #$results{$type}{region_no}=$region_no;
      #$results{$type}{region_total}=$region_total;
      $results{$type}{region_percent} = $pretty_region_percent;
    }
    my $CpG = $results{CpG}{region_percent};
    my $CHG = $results{CHG}{region_percent};
    my $CHH = $results{CHH}{region_percent};
    print
      "$name\t$ref:$start..$end\t$strand\t$strain\t$CpG\t$CHH\t$CHG\t$note\n";

  }
}

