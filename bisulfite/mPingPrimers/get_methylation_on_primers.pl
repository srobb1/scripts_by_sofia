#!/usr/bin/perl -w
use strict;

my $primers = shift;
my $strain = shift; ## A119 NB
open IN, $primers or die "Can't open $primers\n";
my %methy;
while (my $line=<IN>){
  next if $line =~ /ID\sSEQLength/;
  chomp $line;
  #ID      SEQLength       primerSetNum    primerOrient    product_size    start   len     tm      gc%     primerSeq
  #Chr6:6747723..6747922   200     1       left            26      27      68.384  48.148  CGTGCCGTGTATTATTACGCCAATACC
  my ($loc,$seqLen,$primerPairNum,$primerOrient,$productSize,$primerStart,$primerLen,$tm,$gc,$seq) = split /\t/ ,$line;
  next if $productSize !~ /^\d+$/ ;
  next if $primerPairNum > 1 ;
  ##Chr4:11588185..11588384
  my ($ref,$refStart,$refEnd) = $loc =~ /(Chr\d+)\:(\d+)\.\.(\d+)/i;
  my @seq = split '',$seq;
  for (my $i = 0 ; $i < @seq ; $i++){
    my $nt = $seq[$i]; 
    my $primer_pos = $primerStart + $i + 1;
#print "$primerOrient $nt before \n";
    my $ref_pos = $primer_pos  + $refStart;
    next unless $nt eq 'C' and $primerOrient eq 'left' or $nt eq 'G' and $primerOrient eq 'right';
#print $line ," after \n";
    my $type;
    #print "$nt,$ref_pos,$primerOrient\n";;
    foreach my $type (qw(CpG CHH CHG)){
      my ($percent_meth,$total,$meth,$unmeth) = get_percent_methylation($type, $strain ,"$ref:$ref_pos..$ref_pos" ,'/rhome/robb/project/bisulfide/output_all/split_methyl');
#      print "($percent_meth,$total,$meth,$unmeth)\n";
#      if ($percent_meth =~ /^\d/){
      if (defined $meth){
        $methy{$loc}{$ref_pos}{$nt}{$type}{percent}=$percent_meth;
        $methy{$loc}{$ref_pos}{$nt}{$type}{total}=$total;
        $methy{$loc}{$ref_pos}{$nt}{$type}{meth}=$meth;
        $methy{$loc}{$ref_pos}{$nt}{$type}{unmeth}=$unmeth;
        print "$nt $type $strain $ref:$ref_pos..$ref_pos $percent_meth $total $meth $unmeth\n";
      }
      elsif ($percent_meth== -1){
         print $line ," $total \n";
      }
    }
  }
}
sub get_percent_methylation {
  my $type   = shift;    ## CHH, CHG, CpG
  my $strain = shift;
  my $range  = shift;
  my $dir    = shift;##split_methylation

  my ( $ref, $range_start, $range_end ) = $range =~ /(.+)\:(\d+)\.\.(\d+)/;
  my %methy;
  my $file = "$dir/${strain}_${ref}_${type}.split.txt";
  open METHY, $file or return (-1,$file);;    #NB_Chr9_CpG.split.txt
  while ( my $line = <METHY> ) {
    chomp $line;
    my ( $pos, $methyl_yes, $methyl_no ) = split /\t/, $line;
    next unless $pos >= $range_start;
    next unless $pos <= $range_end;
    last if $pos > $range_end;
    $methyl_yes = defined $methyl_yes ? $methyl_yes : 0;
    $methyl_no  = defined $methyl_no  ? $methyl_no  : 0;
    $methy{$pos}{$type}{yes} += $methyl_yes;
    $methy{$pos}{$type}{no}  += $methyl_no;
  }
  my $region_yes = 0;
  my $region_no  = 0;
  foreach my $pos ( sort { $a <=> $b } keys %methy ) {
    foreach my $type ( keys %{ $methy{$pos} } ) {
      my ( $methy_yes, $methy_no ) = ( 0, 0 );
      if ( exists $methy{$pos}{$type}{yes} ) {
        $methy_yes = $methy{$pos}{$type}{yes};
        $region_yes += $methy_yes;
      }
      if ( exists $methy{$pos}{$type}{no} ) {
        $methy_no = $methy{$pos}{$type}{no};
        $region_no += $methy_no;
      }
      if (!exists $methy{$pos}{$type}){
        return ();
      }
      my $total = $methy_yes + $methy_no;
      my $percent = ( $methy_yes / ($total) ) * 100;
  
      #if ( $percent == 0 ) {
      #  $percent = -100;
      #}
      my $pretty_percent = sprintf( '%.1f', $percent );
      return ($pretty_percent ,$total ,$methy_yes ,$methy_no);
    }
  }
  #my $region_total          = $region_yes + $region_no;
  #my $region_percent        = ( $region_yes / $region_total ) * 100;
  #my $pretty_region_percent = sprintf( '%.1f', $region_percent );
  
  #print "$range $pretty_region_percent% $type in $strain\n";
  #print "=========================================================\n";
  #print $to_print;
  #print "=========================================================\n";
  #print "$range $pretty_region_percent% $type in $strain\n";
}
