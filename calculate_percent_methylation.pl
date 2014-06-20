#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#my $dir = shift;
#***A119.A119_p1.fq_bismark_bt2_pe.sam***
#firstBefore, secondIn(notUniq): 387(69),658(94),is:365,ping11(6547):601..5947
#DGGXHXP1:377:C2K44ACXX:8:2212:12680:8801_1:N:0:CCGTCC   99      ping11  387     255     69M     =       658     365     AGGTGTGGAGGTTTGTTGTGTGGAGGTGTTGTGGTTGTTTTTGTTTTTGTCGATTATATTGATATGTGG   =:+ABBBACC;<AF7?AC:AAG;C8?E8CF?D1@FF@@DFFI;FFFII@D@8CEEEAEDD7@BBD;6;:   NM:i:30 XX:Z:3C1C5CCC1C1AC7C1CC1C2CC3CCC1CCCCC1C4C1C1CC2C1C1C2  XM:Z:...z.z.....hxz.h..z.......z.xz.z..xz...hxz.hhhxz.xZ...h.h.xz..h.z.z..      XR:Z:CT XG:Z:CT
#DGGXHXP1:377:C2K44ACXX:8:2212:12680:8801_1:N:0:CCGTCC   147     ping11  658     255     94M     =       387     -365    GACTGAATAAAAAATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGATTTTTATTTTGGTGAAATTCGTCAGCGTTGTTTCCAAGTT  ABAB>A>??ADDDAA>:DDDDD@:>CDEDCCA;6/)CCC=B=B/))*D?*D?90*EEDD9CD>A>EDEE>DCC1AC8A@:IDC8DDA++:=8??  NM:i:7  XX:Z:54G3C2CC8C10C10C   XM:Z:..X.....................H.................................h..hx........h.Z..X..Z..z....HH....h     XR:Z:GA XG:Z:CT

my %code  = qw (X CHG x CHG H CHH h CHH Z CpG z CpG);
my $region = shift ;
my @files = @ARGV; ##goodhits.sam                                     #<$dir/C*txt>;
my $TE    = 'mping';
if ( !@files ) {
  die "include the reads summary from find_good_ping_hits.pl\n";
}
my $i = 0;
foreach my $file (@files) {
  my %methy;
  my %reads;
  my %ref;
  my @path = split /\//, $file;
  my $file_name = pop @path;
  my ( $type, $sample )
    ;    # = $file_name =~ /^(CpG|CHH|CHG)_\w{2,4}_(\w{2,5})\.\2/;
         #print "$type\t$sample\n";
  open IN, $file or die "Can't open $file for reading $!\n";

  while ( my $line = <IN> ) {
    chomp $line;
    ## need to store the lens of all the refs from @SQ
    if ( $line =~ /^\@SQ/ ) {
      #@SQ     SN:mping.Chr9_643091_643093     LN:1636
      my ( $ref, $len ) = $line =~ /SN:(\S+)\s+LN:(\d+)/;
      $ref{$ref} = $len;
      next;
    }
    #@CO sample:RIL60.RIL60_p1.fq_bismark_bt2_pe
    if ( $line =~ /^\@CO sample:(\w+)\.\1_p\d\.fq_bismark/ ) {
      $sample = $1;
      next;
    }
    my @line = split /\t/, $line;
    next unless @line > 10;
    my ( $flag, $ref, $pos, $methy_str ) =
      ( $line[1], $line[2], $line[3], $line[13] );
    $methy_str =~ s/XM:Z://;
    my $strand = ( $flag & 16 ) ? '-' : '+';
    my @methy_str = split '', $methy_str;
    foreach my $methy (@methy_str) {
      my $type = $code{$methy};
      $methy{$sample}{$ref}{$pos}{$methy}++;
      $reads{$sample}{$ref}{$pos}{ $line[0] }++;
      $pos++;
    }
  }

  #print Dumper \%methy;
  my $flank = 600;
  print
    "sample\tref\tregion\tCHH\%\tCpG\%\tCHG\%\tCHH:+;-;T\tCpG:+;-;T\tCHG:+;-;T\tpair-reads\n";
  foreach my $sample ( keys %methy ) {
    foreach my $ref ( keys %{ $methy{$sample} } ) {
      ## 6,547 - 600 = 5947
      ## 6,547 - 603 = 5944
      my $len = $ref{$ref};

      my $end = $len - $flank + 1;
      my %count;
      ## sub region checked 600,300,150
      #for ( my $i = 1 ; $i <= 600 ; $i++ ) { ## for 1-600
      for ( my $i = ($flank-$region+1) ; $i <= $flank ; $i++ ) { ## for 300-600
        if ( exists $methy{$sample}{$ref}{$i} ) {
          foreach my $code ( keys %{ $methy{$sample}{$ref}{$i} } ) {
            my $count = $methy{$sample}{$ref}{$i}{$code};
            foreach my $read ( keys %{ $reads{$sample}{$ref}{$i} } ) {
              $reads{pre}{$sample}{$ref}{$read}++;
            }
            if ( $code ne '.' ) {
              $count{pre}{$code} += $count;
            }
          }
        } else {

          #print "$sample\t$ref\t$i\t.\t0\n";

        }
      }
      #for ( my $i = $end ; $i <= $len ; $i++ ) { ## for flank of 600
      for ( my $i = $end ; $i <= $end+$region+1 ; $i++ ) { ## for flank of 300
        if ( exists $methy{$sample}{$ref}{$i} ) {
          foreach my $code ( keys %{ $methy{$sample}{$ref}{$i} } ) {
            my $count = $methy{$sample}{$ref}{$i}{$code};
            foreach my $read ( keys %{ $reads{$sample}{$ref}{$i} } ) {
              $reads{post}{$sample}{$ref}{$read}++;
            }
            if ( $code ne '.' ) {
              $count{post}{$code} += $count;
            }
          }
        } else {

          #print "$sample\t$ref\t$i\t.\t0\n";

        }
      }
      foreach my $region ( keys %count ) {
        my $reads     = keys %{ $reads{$region}{$sample}{$ref} };
        my $Z         = exists $count{$region}{Z} ? $count{$region}{Z} : 0;
        my $z         = exists $count{$region}{z} ? $count{$region}{z} : 0;
        my $X         = exists $count{$region}{X} ? $count{$region}{X} : 0;
        my $x         = exists $count{$region}{x} ? $count{$region}{x} : 0;
        my $H         = exists $count{$region}{H} ? $count{$region}{H} : 0;
        my $h         = exists $count{$region}{h} ? $count{$region}{h} : 0;
        my $total_H   = $h + $H;
        my $total_X   = $x + $X;
        my $total_Z   = $z + $Z;
        my $percent_H = ($total_H) > 0 ? ( $H / $total_H ) * 100 : 0;
        my $percent_X = ($total_X) > 0 ? ( $X / $total_X ) * 100 : 0;
        my $percent_Z = ($total_Z) > 0 ? ( $Z / $total_Z ) * 100 : 0;
        my $print_H   = sprintf( '%.1f', $percent_H );
        my $print_X   = sprintf( '%.1f', $percent_X );
        my $print_Z   = sprintf( '%.1f', $percent_Z );
        print
          "$sample\t$ref\t$region\t$print_H\t$print_X\t$print_Z\t$H;$h;$total_H\t$X;$x;$total_X\t$Z;$z;$total_Z\t$reads\n";
      }
    }
  }
  #print Dumper \%reads;
}
