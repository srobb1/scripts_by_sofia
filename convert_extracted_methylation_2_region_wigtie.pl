#!/usr/bin/perl -w
use strict;
use Tie::File;
#my $dir = shift;

##head CpG_OB_RIL43.RIL43_p1.fq_bismark_bt2_pe.txt
#Bismark methylation extractor version v0.9.0
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	376063	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	376060	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375985	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375911	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375914	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375985	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1903:2170_1:N:0:CTTGTA	+	Chr1	37611389	Z
#DGGXHXP1:377:C2K44ACXX:8:1101:1903:2170_1:N:0:CTTGTA	+	Chr1	37611384	Z
#DGGXHXP1:377:C2K44ACXX:8:1101:1903:2170_1:N:0:CTTGTA	+	Chr1	37611343	Z

my $file = shift;    #<$dir/C*txt>;
my $range = shift;

use Tie::File;
print "before tie\n";
tie my @array, 'Tie::File', $file or die "...";
print "after tie\n";

my ($range_ref,$range_start,$range_end) = $range =~ /(.+)\:(\d+)\.\.(\d+)/;
my %reads;
  my %methy;
  my @path = split /\//, $file;
  my $file_name = pop @path;
  #my ( $type, $sample ) = $file_name =~ /^(CpG|CHH|CHG)_\w{2,4}_(\w{2,5})\.\2/;
  my ( $type, $sample ) = $file_name =~ /^(CpG|CHH|CHG).*(NB|A119|A123|RIL\d+|HEG4|EG4)/;

  #print "$type\t$sample\n";
#print "before len\n";
 my $array_len = $#array;
  my $i = $array_len /2;
#print "after len\n";
#  my $i=1_000_000; #arbitrary
#while (!defined $array[$i]){
#  $i = int($i)/2;
#}
  my $done = 0;
  my $current_pos;
  my $last_status = 'noneyet';
  my $last_i;
  my $divisor = 2;
  while (!$done){
    print "step1: $i\n";
    my $status;
    my $line = $array[int($i)];
 print "step2: at arr[$i] is $line\n";
    my ( $read, $strand, $ref, $pos, $methy ) = split /\t/, $line;
    $current_pos = $pos;
    $done = ((($range_start-$pos)> 0) and (($range_start-$pos)<1000)) ? 1 : 0; 
    last if $done;
    if ($pos < $range_start ){
      $status = 'too_small';
      $i = $i + $i/$divisor;

    print "step3: ",int($i),"\n";
    }elsif ($pos > $range_start){
      $status = 'too_big';
      $i = $i - ($i/$divisor);
    print "step3: ",int($i),"\n";
    }else {
print "noneof the above\n";
     } 
    print "$range_start..$range_end array[",int$i,"] is $pos $status\n";
    if ($last_status eq 'too_small' and $status eq 'too_big'){
      $divisor = $divisor * 2;
      ## reset to last $i and last status
      $i = $last_i;
$i = $i - ($i/$divisor);
      $status = $last_status;
      print "reset last_i=$i last_status=$status\n";
    } 
    $last_status = $status;
    $last_i = $i;
  }
  while ($current_pos < $range_end +1){
    my $line = $array[$i];
    my ( $read, $strand, $ref, $pos, $methy ) = split /\t/, $line;
    $current_pos = $pos;
    $strand = 0;
    $methy{$sample}{$ref}{$pos}{$strand}{$type}{$methy}++;
    $i++;
  }

  foreach my $sample ( keys %methy ) {
    foreach my $ref ( keys %{ $methy{$sample} } ) {
      foreach my $pos ( sort { $a <=> $b } keys %{ $methy{$sample}{$ref} } ) {
        $i++;
        foreach my $strand ( keys %{ $methy{$sample}{$ref}{$pos} } ) {
          foreach my $type ( keys %{ $methy{$sample}{$ref}{$pos}{$strand} } ) {
            my @methy_codes =
              keys %{ $methy{$sample}{$ref}{$pos}{$strand}{$type} };
            my $code   = $methy_codes[0];
            my $big    = uc $code;
            my $little = lc $code;
            my ( $methy_yes, $methy_no ) = ( 0, 0 );
            if ( exists $methy{$sample}{$ref}{$pos}{$strand}{$type}{$big} ) {
              $methy_yes = $methy{$sample}{$ref}{$pos}{$strand}{$type}{$big};
            }
            if ( exists $methy{$sample}{$ref}{$pos}{$strand}{$type}{$little} ) {
              $methy_no = $methy{$sample}{$ref}{$pos}{$strand}{$type}{$little};
            } 
              my $total          = $methy_yes + $methy_no;
              my $percent        = ( $methy_yes / ($total) ) * 100;
              if ($percent == 0){
                $percent = -100;
              }
              my $pretty_percent = sprintf( '%.1f', $percent );
              print
                "$type\t$ref\t",$pos-1,"\t$pos\t$pretty_percent\t$total\t$methy_yes\t$methy_no\n";
          }
        }
      }
    }
  }
