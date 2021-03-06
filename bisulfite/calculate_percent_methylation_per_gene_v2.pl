#!/usr/bin/perl -w
use strict;
use Bio::DB::SeqFeature::Store;
use Data::Dumper;

my $sqlite = shift;
my $dir    = shift ; #"/rhome/robb/project/bisulfide/output_all/split_methyl";    #shift
my @strains = @ARGV;

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
foreach my $strain (@strains ) {
#foreach my $strain (qw (A119 NB RIL16 RIL25 RIL43 RIL60)) {
  foreach my $ref (
    qw (Chr1 Chr2 Chr3 Chr4 Chr5 Chr6 Chr7 Chr8 Chr9 Chr10 Chr11 Chr12 ChrUn ChrSy)
    )
  {
    my %results;
    my @features = $db_obj->features( -seq_id => $ref, -types => ['gene'] );
    foreach my $type (qw (CHH CpG CHG)) {
      my %methy;
      open IN, "$dir/${strain}_${ref}_${type}.split.txt"
        or die "Can't open $dir/${strain}_${ref}_${type}.split.txt";
      while ( my $line = <IN> ) {
        chomp $line;
        my ( $pos, $methyl_yes, $methyl_no ) = split /\t/, $line;
        $methyl_yes = defined $methyl_yes ? $methyl_yes : 0;
        $methyl_no  = defined $methyl_no  ? $methyl_no  : 0;
        $methy{$pos}{yes} += $methyl_yes;
        $methy{$pos}{no}  += $methyl_no;
      }
      #my @features = $db_obj->get_features_by_location( $ref );
      #my @features = $db_obj->get_features_by_type('gene');
      foreach my $f (@features) {
        my $start      = $f->start;
        my $end        = $f->end;
        my $strand     = $f->strand;
        my %attr       = $f->attributes;
        my $note       = ${ $attr{'Note'} }[0];
        my $name       = $f->name;
        my $region_yes = 0;
        my $region_no  = 0;
#        print join( '--', $name, $f->ref, $start, $end, $strand, $note ), "\n";

        for ( my $pos = $start ; $pos <= $end ; $pos++ ) {
          next unless exists $methy{$pos};
          my ( $methy_yes, $methy_no ) = ( 0, 0 );
          if ( exists $methy{$pos}{yes} ) {
            $methy_yes = $methy{$pos}{yes};
            $region_yes += $methy_yes;
          }
          if ( exists $methy{$pos}{no} ) {
            $methy_no = $methy{$pos}{no};
            $region_no += $methy_no;
          }

        }
        my $region_total = $region_yes + $region_no;
        my $region_percent;
        my $pretty_region_percent;
        if ( $region_total == 0 or !defined $region_total ) {
          $pretty_region_percent = 'n/a';
        } else {
          $region_percent = ( $region_yes / $region_total ) * 100;
          $pretty_region_percent = sprintf( '%.1f', $region_percent );
        }

        $results{$name}{region}                = "$ref:$start..$end";
        $results{$name}{region_percent}{$type} = $pretty_region_percent;
        $results{$name}{strand}                = $strand;
        $results{$name}{note}                  = $note;
      }
    }
    foreach my $name ( keys %results ) {
      my $CpG    = $results{$name}{region_percent}{CpG};
      my $CHG    = $results{$name}{region_percent}{CHG};
      my $CHH    = $results{$name}{region_percent}{CHH};
      my $region = $results{$name}{region};
      my $strand = $results{$name}{strand};
      my $note   = $results{$name}{note};
      print "$name\t$region\t$strand\t$strain\t$CpG\t$CHH\t$CHG\t$note\n";
    }
  }
}
