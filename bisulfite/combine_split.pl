#!/usr/bin/perl -w
use strict;

my $dir = shift;
my @files = <$dir/*split.txt>;
my %meth;
my %strains;
my %locs;
foreach my $file (@files){
  #A119_2_mping.Chr4_27129772_27129774_CHG.split.txt
  my ($strain,$loc,$type) = $file =~ /(\w+_\d)_(.+)_(C..)\.split\.txt/;
  #print "($strain,$loc,$type)\n";
  $strains{$strain}++;
  open IN, $file or die "can't open $file\n";
  while (my $line = <IN>){
    my ($pos,$meth_count,$no_meth_count) = split /\s+/ ,$line;
    #print "\t($pos,$meth_count,$no_meth_count)\n";
    $meth{$loc}{$pos}{$type}{$strain}{meth} =$meth_count;
    $meth{$loc}{$pos}{$type}{$strain}{nometh} =$no_meth_count;
  }
}
my @header =('loc','insert','pos');
foreach my $strain (sort keys %strains){
  push @header, "$strain:CHH","$strain:chh","$strain:%CHH","$strain:CHG","$strain:chg","$strain:%CHG","$strain:CpG","$strain:cpg","$strain:%CpG";
} 
print join ("\t",@header), "\n";

foreach my $loc (sort keys %meth){
  my $general_loc = $loc;
  foreach my $pos (sort {$a<=>$b} keys %{$meth{$loc}}){
    my $empty = 'insert';
    if ($loc =~ /empty/){
      $empty = 'empty_site';
      $general_loc =~ s/empty\.//;
    } 
    my @row = ($general_loc,$empty,$pos);
    my ($CHH,$chh,$CHH_per,  $CHG,$chg,$CHG_per, $CpG,$cpg,$CpG_per) = ('-','-','-','-','-','-','-','-','-');
    #my ($CHH,$chh,$CHH_per,  $CHG,$chg,$CHG_per, $CpG,$cpg,$CpG_per) = ('','','','','','','','','');
    foreach my $strain (sort keys %strains){
      if (exists $meth{$loc}{$pos}{CHH}{$strain}){
        $CHH =  $meth{$loc}{$pos}{CHH}{$strain}{meth}; 
        $chh =  $meth{$loc}{$pos}{CHH}{$strain}{nometh}; 
        $CHH_per =  $CHH / ($CHH + $chh);
      }
      if (exists $meth{$loc}{$pos}{CHG}{$strain}){
        $CHG = $meth{$loc}{$pos}{CHG}{$strain}{meth}; 
        $chg = $meth{$loc}{$pos}{CHG}{$strain}{nometh}; 
        $CHG_per =  $CHG / ($CHG + $chg);
      }
      if (exists $meth{$loc}{$pos}{CpG}{$strain}){
        $CpG = $meth{$loc}{$pos}{CpG}{$strain}{meth}; 
        $cpg = $meth{$loc}{$pos}{CpG}{$strain}{nometh}; 
        $CpG_per =  $CpG / ($CpG + $cpg);
      }
    
      push @row,  $CHH, $chh, $CHH_per, $CHG, $chg, $CHG_per, $CpG, $cpg, $CpG_per;
    }   
    print join ("\t",@row), "\n";
  }
}

