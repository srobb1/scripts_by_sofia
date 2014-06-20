#!/usr/bin/perl -w
use strict;
use DBI;
my $db_file    = 'bismark';
my $range_type = shift;       ## CHH, CHG, CpG
my $strain     = shift;
my $range      = shift;

my $dbh = DBI->connect( "dbi:SQLite:dbname=$db_file", "", "" );
my ( $range_ref, $range_start, $range_end ) = $range =~ /(.+)\:(\d+)\.\.(\d+)/;
my %reads;
my %methy;

my $sql =
  qq (SELECT s.name,t.name,m.ref,m.pos,m.yes, m.no from METHYLATION m, TYPES t, STRAINS s WHERE m.ref = "$range_ref" and m.pos >= $range_start and m.pos <= $range_end and t.name = "$range_type" and t.ID = m.typeID and s.name="$strain" and s.id = m.strainID);

print $sql,"\n";
__END__
my $sth = $dbh->prepare($sql);
$sth->execute;

while ( my $ary_ref = $sth->fetchrow_arrayref ) {
  my ( $sample, $type, $ref, $pos, $methyl_yes, $methyl_no ) = @$ary_ref;
  $methyl_yes = defined $methyl_yes ? $methyl_yes : 0;
  $methyl_no = defined $methyl_no ? $methyl_no : 0;
  $methy{$sample}{$ref}{$pos}{$type}{$methyl_yes}++;
  $methy{$sample}{$ref}{$pos}{$type}{$methyl_no}++;
}
foreach my $sample ( keys %methy ) {
  foreach my $ref ( keys %{ $methy{$sample} } ) {
    foreach my $pos ( sort { $a <=> $b } keys %{ $methy{$sample}{$ref} } ) {
      foreach my $type ( keys %{ $methy{$sample}{$ref}{$pos} } ) {
        my ( $methy_yes, $methy_no ) = ( 0, 0 );
        if ( exists $methy{$sample}{$ref}{$pos}{$type}{1} ) {
          $methy_yes = $methy{$sample}{$ref}{$pos}{$type}{1};
        }
        if ( exists $methy{$sample}{$ref}{$pos}{$type}{0} ) {
          $methy_no = $methy{$sample}{$ref}{$pos}{$type}{0};
        }
        my $total = $methy_yes + $methy_no;
        my $percent = ( $methy_yes / ($total) ) * 100;
        if ( $percent == 0 ) {
          $percent = -100;
        }
        my $pretty_percent = sprintf( '%.1f', $percent );
        print "$type\t$ref\t", $pos - 1,
          "\t$pos\t$pretty_percent\t$total\t$methy_yes\t$methy_no\n";
      }
    }
  }
}
