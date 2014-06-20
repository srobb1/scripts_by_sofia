#!/usr/bin/perl -w
use strict;
use Data::Dumper;
#my $dir = shift;
#***A119.A119_p1.fq_bismark_bt2_pe.sam***
#firstBefore, secondIn(notUniq): 387(69),658(94),is:365,ping11(6547):601..5947
#DGGXHXP1:377:C2K44ACXX:8:2212:12680:8801_1:N:0:CCGTCC   99      ping11  387     255     69M     =       658     365     AGGTGTGGAGGTTTGTTGTGTGGAGGTGTTGTGGTTGTTTTTGTTTTTGTCGATTATATTGATATGTGG   =:+ABBBACC;<AF7?AC:AAG;C8?E8CF?D1@FF@@DFFI;FFFII@D@8CEEEAEDD7@BBD;6;:   NM:i:30 XX:Z:3C1C5CCC1C1AC7C1CC1C2CC3CCC1CCCCC1C4C1C1CC2C1C1C2  XM:Z:...z.z.....hxz.h..z.......z.xz.z..xz...hxz.hhhxz.xZ...h.h.xz..h.z.z..      XR:Z:CT XG:Z:CT
#DGGXHXP1:377:C2K44ACXX:8:2212:12680:8801_1:N:0:CCGTCC   147     ping11  658     255     94M     =       387     -365    GACTGAATAAAAAATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGATTTTTATTTTGGTGAAATTCGTCAGCGTTGTTTCCAAGTT  ABAB>A>??ADDDAA>:DDDDD@:>CDEDCCA;6/)CCC=B=B/))*D?*D?90*EEDD9CD>A>EDEE>DCC1AC8A@:IDC8DDA++:=8??  NM:i:7  XX:Z:54G3C2CC8C10C10C   XM:Z:..X.....................H.................................h..hx........h.Z..X..Z..z....HH....h     XR:Z:GA XG:Z:CT

my %code = qw (X CHG x CHG H CHH h CHH Z CpG z CpG);
my $file =  shift; #<$dir/C*txt>;
#my $range = shift;
#my ($range_ref,$range_start,$range_end) = $range =~ /(.+)\:(\d+)\.\.(\d+)/;
my $i=0;
  my %methy;
  my @path = split /\// , $file;
  my $file_name = pop @path; 
  #RIL60.RIL60_p1.fq_bismark_bt2_pe.sam
  my ($sample) = $file_name =~ /(\w+).+sam/;
  open IN , $file or die "Can't open $file for reading $!\n";
  while (my $line = <IN>){
    next if $line =~ /^\@SQ/;
    chomp $line;
    my @line = split /\t/ , $line;
    next unless @line > 10;  
    my ($flag,$ref,$pos,$methy_str) = ($line[1],$line[2],$line[3],$line[13]); 
    #next unless $ref eq $range_ref;
    #next if $pos > $range_end;
    $methy_str =~ s/XM:Z://;
    #my $len = length $methy_str;
    #next if $pos+$len < $range_start;
    my $strand = ($flag & 16) ? '-' : '+';
    my @methy_str = split '',$methy_str;
    $strand = 0;
    foreach my $methy (@methy_str){
      if ($methy ne '.'){
        my $type = $code{$methy};
        $methy{$sample}{$type}{$ref}{$pos}{$strand}{$methy}++ if defined $type;
      }
      $pos++;
    }
    
  }
  foreach my $sample (keys %methy){
  foreach my $type (keys %{$methy{$sample}}){
    open OUTWIG ,">$type"."_$sample.justPings.wig" or die "Can't open  $sample.$type.wig\n";
    print OUTWIG "track type=wiggle_0 name=\"$sample.$type\" description=\"$type in $sample\"\n";
    foreach my $ref (keys %{$methy{$sample}{$type}}){
    foreach my $pos (sort {$a <=> $b} keys %{$methy{$sample}{$type}{$ref}}){
      $i++;
      foreach my $strand (keys  %{$methy{$sample}{$type}{$ref}{$pos}}){
          my @methy_codes = keys  %{$methy{$sample}{$type}{$ref}{$pos}{$strand}};
          my $code = $methy_codes[0];
          my $big = uc $code;
          my $little = lc $code;
          my ($methy_yes,$methy_no)=(0,0);
          if (exists $methy{$sample}{$type}{$ref}{$pos}{$strand}{$big}){
            $methy_yes = $methy{$sample}{$type}{$ref}{$pos}{$strand}{$big}; 
          }if(exists $methy{$sample}{$type}{$ref}{$pos}{$strand}{$little}){
            $methy_no = $methy{$sample}{$type}{$ref}{$pos}{$strand}{$little};
          }
            my $total = $methy_yes+$methy_no;
            my $percent = ($methy_yes/($total)) *100;
            my $pretty_percent   = sprintf( '%.1f', $percent );
            my $wig_s = $pos-1;
            #wig format
            #ping.Chr1_33282257_33282259     237     238     0.0
            #ping.Chr1_33282257_33282259     257     258     0.0
            print OUTWIG "$ref\t$wig_s\t$pos\t$pretty_percent;$total;$methy_yes;$methy_no\n";
        }
      } 
    }
  }
  }
