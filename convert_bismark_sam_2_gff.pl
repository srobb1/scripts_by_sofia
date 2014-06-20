#!/usr/bin/perl -w
use strict;
use Data::Dumper;
#my $dir = shift;
#***A119.A119_p1.fq_bismark_bt2_pe.sam***
#firstBefore, secondIn(notUniq): 387(69),658(94),is:365,ping11(6547):601..5947
#DGGXHXP1:377:C2K44ACXX:8:2212:12680:8801_1:N:0:CCGTCC   99      ping11  387     255     69M     =       658     365     AGGTGTGGAGGTTTGTTGTGTGGAGGTGTTGTGGTTGTTTTTGTTTTTGTCGATTATATTGATATGTGG   =:+ABBBACC;<AF7?AC:AAG;C8?E8CF?D1@FF@@DFFI;FFFII@D@8CEEEAEDD7@BBD;6;:   NM:i:30 XX:Z:3C1C5CCC1C1AC7C1CC1C2CC3CCC1CCCCC1C4C1C1CC2C1C1C2  XM:Z:...z.z.....hxz.h..z.......z.xz.z..xz...hxz.hhhxz.xZ...h.h.xz..h.z.z..      XR:Z:CT XG:Z:CT
#DGGXHXP1:377:C2K44ACXX:8:2212:12680:8801_1:N:0:CCGTCC   147     ping11  658     255     94M     =       387     -365    GACTGAATAAAAAATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGATTTTTATTTTGGTGAAATTCGTCAGCGTTGTTTCCAAGTT  ABAB>A>??ADDDAA>:DDDDD@:>CDEDCCA;6/)CCC=B=B/))*D?*D?90*EEDD9CD>A>EDEE>DCC1AC8A@:IDC8DDA++:=8??  NM:i:7  XX:Z:54G3C2CC8C10C10C   XM:Z:..X.....................H.................................h..hx........h.Z..X..Z..z....HH....h     XR:Z:GA XG:Z:CT

my %code = qw (X CHG x CHG H CHH h CHH Z CpG z CpG);
my @files = @ARGV; #<$dir/C*txt>;
my $i=0;
foreach my $file (@files){
  my %methy;
  my @path = split /\// , $file;
  my $file_name = pop @path; 
  my ($type,$sample);# = $file_name =~ /^(CpG|CHH|CHG)_\w{2,4}_(\w{2,5})\.\2/;
  #print "$type\t$sample\n";
  open IN , $file or die "Can't open $file for reading $!\n";
  while (my $line = <IN>){
    next if $line =~ /^\@SQ/;
    chomp $line;
    #if ($line =~ /\*\*\*(\w+)\.\1_p\d\.fq_bismark_bt2_pe.sam/){
    if ($line =~ /^\@CO sample:(\w+)\.\1_p\d\.fq_bismark/ ) {
      $sample = $1;
      next;
    }
    my @line = split /\t/ , $line;
    next unless @line > 10;  
    my ($flag,$ref,$pos,$methy_str) = ($line[1],$line[2],$line[3],$line[13]); 
    $methy_str =~ s/XM:Z://;
    my $strand = ($flag & 16) ? '-' : '+';
    my @methy_str = split '',$methy_str;
    $strand = 0;
    foreach my $methy (@methy_str){
      if ($methy ne '.'){
        my $type = $code{$methy};
        $methy{$sample}{$ref}{$pos}{$strand}{$type}{$methy}++ if defined $type;
      }
      $pos++;
    }
    
  }
  foreach my $sample (keys %methy){
  foreach my $ref (keys %{$methy{$sample}}){
    foreach my $pos (sort {$a <=> $b} keys %{$methy{$sample}{$ref}}){
      $i++;
      foreach my $strand (keys  %{$methy{$sample}{$ref}{$pos}}){
        foreach my $type (keys %{$methy{$sample}{$ref}{$pos}{$strand}} ){
          my @methy_codes = keys  %{$methy{$sample}{$ref}{$pos}{$strand}{$type}};
          if (scalar @methy_codes == 1 and lc($methy_codes[0]) eq  $methy_codes[0]){
            ## if only one methy (z|Z) is reported, it must be upppercased version for it to be printed, otherwise print both
            next;
          }
          my $code = $methy_codes[0];
          my $big = uc $code;
          my $little = lc $code;
          my ($methy_yes,$methy_no)=(0,0);
          if (exists $methy{$sample}{$ref}{$pos}{$strand}{$type}{$big}){
            $methy_yes = $methy{$sample}{$ref}{$pos}{$strand}{$type}{$big}; 
          }if(exists $methy{$sample}{$ref}{$pos}{$strand}{$type}{$little}){
            $methy_no = $methy{$sample}{$ref}{$pos}{$strand}{$type}{$little};
          }
          if ($methy_yes){
            my $total = $methy_yes+$methy_no;
            my $percent = ($methy_yes/($total)) *100;
            my $pretty_percent   = sprintf( '%.1f', $percent );
            print "$ref\t$sample.$type\tmethylated_cytosine\t$pos\t$pos\t$pretty_percent\t$strand\t.\tID=$sample.$ref$strand.$pos.$type.$i;Name=$ref.$pos.$type;Note=$type $big:$methy_yes $little:$methy_no total:$total percent:$pretty_percent;sample=$sample;\n";
          }
        }
      } 
    }
  }
  }
} 
