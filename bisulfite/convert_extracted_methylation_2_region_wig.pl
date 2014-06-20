#!/usr/bin/perl -w
use strict;

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
my ($range_ref,$range_start,$range_end) = $range =~ /(.+)\:(\d+)\.\.(\d+)/;
my %reads;
my $i     = 0;
  my %methy;
  my @path = split /\//, $file;
  my $file_name = pop @path;
  #my ( $type, $sample ) = $file_name =~ /^(CpG|CHH|CHG)_\w{2,4}_(\w{2,5})\.\2/;
  my ( $type, $sample ) = $file_name =~ /^(CpG|CHH|CHG).*(NB|A119|A123|RIL\d+|HEG4|EG4)/;

  #print "$type\t$sample\n";
  open IN, $file or die "Can't open $file for reading $!\n";
  while ( my $line = <IN> ) {
    chomp $line;
    next if $line =~ /^Bismark methylation extractor version/;
    my ( $read, $strand, $ref, $pos, $methy ) = split /\t/, $line;
    next unless $ref eq $range_ref;
    next unless $range_start <= $pos; 
    next unless $range_end >= $pos; 
    last if $pos > $range_end; 

    #next unless uc($methy) eq $methy;
    $strand = 0;
    $methy{$sample}{$ref}{$pos}{$strand}{$type}{$methy}++;
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
