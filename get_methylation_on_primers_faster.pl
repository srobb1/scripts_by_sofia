#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $primers = shift;
my $strain  = shift;    ## A119 NB
open IN, $primers or die "Can't open $primers\n";
my %methy;
while ( my $line = <IN> ) {
  next if $line =~ /ID\sSEQLength/;
  chomp $line;

#ID      SEQLength       primerSetNum    primerOrient    product_size    start   len     tm      gc%     primerSeq
#Chr6:6747723..6747922   200     1       left            26      27      68.384  48.148  CGTGCCGTGTATTATTACGCCAATACC
  my (
       $loc,         $seqLen,      $primerPairNum, $primerOrient,
       $productSize, $primerStart, $primerLen,     $tm,
       $gc,          $seq
  ) = split /\t/, $line;
  next if $productSize !~ /^\d+$/;
  next if $primerPairNum > 1;
  ##Chr4:11588185..11588384
  my ( $ref, $refStart, $refEnd ) = $loc =~ /^(Chr\d+)\:(\d+)\.\.(\d+)/i;
  my $primer_ref_start   = $primerStart + $refStart;
  my $primer_ref_end   =   $primer_ref_start + $primerLen - 1;
  my @seq = split '', $seq;
  foreach my $type (qw(CpG CHH CHG)) {
    my $calculated = get_percent_methylation( $type, $strain, "$ref:$primer_ref_start..$primer_ref_end",
                      '/rhome/robb/project/bisulfide/output_all/split_methyl' );
    foreach my $ref_pos ( sort { $a <=> $b } keys %$calculated ) {
      my $percent_meth = ${$calculated}{$ref_pos}{pretty_percent};
      my $total        = ${$calculated}{$ref_pos}{total};
      my $meth         = ${$calculated}{$ref_pos}{methy_yes};
      my $unmeth       = ${$calculated}{$ref_pos}{methy_no};
      my $i            = $ref_pos - $primer_ref_start;
      my $nt           = $seq[$i];
#print "$primerOrient $nt\n";
      next
        unless $nt eq 'C' and $primerOrient eq 'left'
          or $nt eq 'G' and $primerOrient eq 'right';
      #my $type;

      #print "$nt,$ref_pos,$primerOrient\n";;
      if ( defined $meth ) {
        $methy{$loc}{$ref_pos}{$nt}{$type}{percent} = $percent_meth;
        $methy{$loc}{$ref_pos}{$nt}{$type}{total}   = $total;
        $methy{$loc}{$ref_pos}{$nt}{$type}{meth}    = $meth;
        $methy{$loc}{$ref_pos}{$nt}{$type}{unmeth}  = $unmeth;
        print
          "$nt $type $strain $ref:$ref_pos..$ref_pos $percent_meth $total $meth $unmeth\n";
      } elsif ( $percent_meth == -1 ) {
        print $line , " $total \n";
      }
    }
  }
}

sub get_percent_methylation {
  my $type   = shift;    ## CHH, CHG, CpG
  my $strain = shift;
  my $range  = shift;
  my $dir    = shift;    ##split_methylation

  my ( $ref, $range_start, $range_end ) = $range =~ /(.+)\:(\d+)\.\.(\d+)/;
  my %methy;
  my $file = "$dir/${strain}_${ref}_${type}.split.txt";
  open METHY, $file or warn "Can't open $file\n";
  #A119_Chr7:22059928..22060127|mPing@A119_2.Chr7_CHH.split.txt
  while ( my $line = <METHY> ) {
    chomp $line;
    my ( $pos, $methyl_yes, $methyl_no ) = split /\t/, $line;
    next
      unless $pos >= $range_start;
    next
      unless $pos <= $range_end;
    last
      if $pos > $range_end;
    $methyl_yes = defined $methyl_yes ? $methyl_yes : 0;
    $methyl_no  = defined $methyl_no  ? $methyl_no  : 0;
    $methy{$pos}{$type}{yes} += $methyl_yes;
    $methy{$pos}{$type}{no}  += $methyl_no;
  }
  my $region_yes = 0;
  my $region_no  = 0;
  my %calculated;
  foreach my $pos ( sort { $a <=> $b } keys %methy ) {
    foreach my $type ( keys %{ $methy{$pos} } ) {
      my ( $methy_yes, $methy_no ) = ( 0, 0 );
      if ( exists $methy{$pos}{$type} ) {
        if ( exists $methy{$pos}{$type}{yes} ) {
          $methy_yes = $methy{$pos}{$type}{yes};
        }
        if ( exists $methy{$pos}{$type}{no} ) {
          $methy_no = $methy{$pos}{$type}{no};
        }
        my $total = $methy_yes + $methy_no;
        $calculated{$pos}{total} = $total;
        my $percent = ( $methy_yes / ($total) ) * 100;
        $calculated{$pos}{pretty_percent} = sprintf( '%.1f', $percent );
        $calculated{$pos}{methy_yes}      = $methy_yes;
        $calculated{$pos}{methy_no}       = $methy_yes;
      } elsif ( !exists $methy{$pos}{$type} ) {
        $calculated{$pos}{methy_yes} = '';
        $calculated{$pos}{methy_no}  = '';
      }
    }
  }
  return ( \%calculated );
}
