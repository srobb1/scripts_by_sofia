#!/usr/bin/perl -w
use strict;
use Data::Dumper;
## run 'cat file.vcf | vcf-to-tab > out.tab' before running this script

my $vcf_tab = shift;
open VCFTAB, "$vcf_tab" or die;


##CHROM  POS     REF     EG4_2   HEG4_2
#Chr1    1117    A       A/A     A/C
chomp (my $header = <VCFTAB>);
my @header = split /\t/ , $header;
shift @header;
shift @header;
shift @header;
my @strains = @header;
my $strain_count = @strains;
my %SNPs;
my %ref;
my %binary;
while (my $line = <VCFTAB>){
  chomp $line;
  ##if you do not want to use any line without all strains having full data use next line
  next if $line =~ /\.\/\./; #skip whole line if ./. 
  next if $line =~/threw out/;
  my ($chr , $pos , $ref_nt , @snp) = split /\t/ , $line;
  for (my $i=0 ; $i < $strain_count ; $i++){
    $SNPs{"$chr.$pos"}{$header[$i]}=$snp[$i];
    $ref{"$chr.$pos"}=$ref_nt;
    for (my $i=0 ; $i < @snp ; $i++){
      my ($a1 , $a2) = $snp[$i] =~ /(.)\/(.)/;
      #every het position wil be diff from ref
      if ($a1 ne $a2){
         $binary{"$chr.$pos"}{$strains[$i]} = 2; ##hets ok, note as 2
      }
      #not enough data
      elsif ($a1 eq '.'){ ##no 9 , count as same as ref
         $binary{"$chr.$pos"}{$strains[$i]} = 9;
      }
      #homozygous position that is diff from ref
      elsif ( ($a1 eq $a2) and ($ref_nt ne $a2)){
         $binary{"$chr.$pos"}{$strains[$i]} = 1;
      }
      #homozygouos position that is same as ref 
      elsif (($a1 eq $a2) and ($ref_nt eq $a2)) {
         $binary{"$chr.$pos"}{$strains[$i]} = 0;
     }
    }
  }
}

my %counts;
my %hets;
my %diff_from_ref;
my %diff_bwt_indiv;
foreach my $loc (sort keys %SNPs){
  #foreach my $pos (sort {$a <=> $b} keys %{$SNPs{$chr}}){
      my @indiv = sort keys %{$SNPs{$loc}};
      for (my $i=0 ; $i < @indiv ; $i++){
        my $indiv_1 = $indiv[$i];
        my $snp = $SNPs{$loc}{$indiv_1};
        my ($a1 , $a2) = $snp =~ /(.)\/(.)/;
        #if ($a1 ne $2){
        #  $hets{$indiv_1}++;
        #  #next; ## no hets
        #}
        next if $a1 eq '.';
        ## if alleles are different, count it
        my ($chr,$pos) = $loc =~ /(\w+)\.(\d+)/;
        if ($a1 ne $a2){
          $diff_from_ref{$indiv_1}{het}{count}++;
          $diff_from_ref{$indiv_1}{het}{loc}{$loc} = $snp;
        }
        ##these will not be hets
        elsif ($a1 ne $ref{"$chr.$pos"} ){
          $diff_from_ref{$indiv_1}{diff}{count}++;
          $diff_from_ref{$indiv_1}{diff}{loc}{$loc} = $snp
        }else{
          $diff_from_ref{$indiv_1}{same}{count}++;
          $diff_from_ref{$indiv_1}{same}{loc}{$loc} = $snp;
        }

        for (my $j=0 ; $j < @indiv ; $j++){ 
           my $indiv_2 = $indiv[$j];
           next if $indiv_1 eq $indiv_2;
           next if exists $counts{$indiv_2}{$indiv_1};
           next if $SNPs{$loc}{$indiv_2} eq './.';
           if ( $SNPs{$loc}{$indiv_1} ne $SNPs{$loc}{$indiv_2}){
             $diff_bwt_indiv{$loc}{"$indiv_1,$indiv_2"}= "$SNPs{$loc}{$indiv_1},$SNPs{$loc}{$indiv_2}";
             $counts{$indiv_1}{$indiv_2}{diff}++;
           } else {
             $counts{$indiv_1}{$indiv_2}{same}++;
           }
        }           
      }
  #}
}
print "strain comparisons\n";
print "comparison\tSNP differences between strain and ref\n";
foreach my $i_1 ( sort keys %counts){
  print "ref/$i_1\t$diff_from_ref{$i_1}{diff}{count}\n";
}
print "\n\ncomparison\tNT the same as ref\n";
foreach my $i_1 ( sort keys %counts){
  print "ref/$i_1\t$diff_from_ref{$i_1}{same}{count}\n";
}
print "\n\ncomparison\tSNP heterozygous: these snps not included in above counts\n";
foreach my $i_1 ( sort keys %counts){
  print $i_1,"(hets)\t$diff_from_ref{$i_1}{het}{count}\n";
}

print "\n\ncomparison\tSNP differences between strains (hets included) 'A/A'='A/A', 'A/T' = 'A/T', 'A/A' ne 'A/T'\n";
foreach my $i_1 ( sort keys %counts){
  foreach my $i_2 (sort keys %{$counts{$i_1}}){
    print "$i_1/$i_2\t$counts{$i_1}{$i_2}{diff}\n";
  }  
}
print "\n\ncomparison\tSNP shared (hets included) 'A/A' = 'A/A' , 'A/T' = 'A/T', 'A/A' ne 'A/T'\n";
foreach my $i_1 ( sort keys %counts){
  foreach my $i_2 (sort keys %{$counts{$i_1}}){
    print "$i_1/$i_2\t$counts{$i_1}{$i_2}{same}\n";
  }
}  

my %shared;
foreach my $loc (sort keys %binary){
  my %haveSNP;
  foreach my $strain (keys %{$binary{$loc}}){
    if ($binary{$loc}{$strain} == 1){
      push @{$haveSNP{homo}} , $strain;
    }elsif($binary{$loc}{$strain} == 2){
      push @{$haveSNP{het}} , $strain;
    }
  }
  if (exists $haveSNP{homo}){
    my $haveSNP_homo = join (',', sort @{$haveSNP{homo}});
    $shared{$haveSNP_homo}{homo}++;
  }
  if (exists $haveSNP{het}){
    my $haveSNP_het = join (',', sort @{$haveSNP{het}});
    $shared{$haveSNP_het}{het}++;
  }
}

print "\n\nnumber of Homozygous SNPs that are unique to this set of strains\n";
foreach my $strains (sort keys %shared){
  my $count = $shared{$strains}{homo};
  print "$strains\t$count\n" if defined $count;
}

print "\n\nnumber of Heterozygous SNPs that are unique in this set of strains\n";
foreach my $strains(sort  keys %shared){
  my $count = $shared{$strains}{het};
  print "$strains\t$count\n" if defined $count;
}
## this is a strict comparison of the strains, no comparison to REF, therefore some strains may be same
## as ref and not a SNP
my %shared_NT;
foreach my $loc (sort keys %SNPs){
  my %haveSNP;
  foreach my $strain (sort keys %{$SNPs{$loc}}){
    my $snp = $SNPs{$loc}{$strain};
    $haveSNP{$snp}{$strain}=1;
  }
  foreach my $snp (sort keys %haveSNP){
     my $strains = join (',' , sort (keys %{$haveSNP{$snp}}));
     my ($a1,$a2) = split /\// , $snp;
     #print "$strains--$a1--$a2\n";
     my $type = $a1 eq $a2 ? 'homo' : 'het' ;
     $shared_NT{$type}{$strains}++;
  }
}

#print "\n\nnumber of Homozygous NT that are unique to this set of strains\n";
#foreach my $strains (sort keys %{$shared_NT{homo}}){
#  my $count = $shared_NT{homo}{$strains};
#  print "$strains\t$count\n" if defined $count;
#}

#print "\n\nnumber of Heterozygous NT that are unique in this set of strains\n";
#foreach my $strains(sort  keys %{$shared_NT{het}}){
#  my $count = $shared_NT{het}{$strains};
#  print "$strains\t$count\n" if defined $count;
#}




print "\n\n***location of homozygous SNPs\n";
foreach my $strain(sort keys %diff_from_ref){
  foreach my $loc ( keys %{$diff_from_ref{$strain}{diff}{loc}}){
    print "$strain\t$loc\n";
  }
}

print "\n\n***location of heterozygous SNPs\n";
foreach my $strain(sort keys %diff_from_ref){
  foreach my $loc ( keys %{$diff_from_ref{$strain}{het}{loc}}){
    print "$strain\t$loc\n";
  }
}

##  $diff_bwt_indiv{$loc}{"$indiv_1,$individ_2"}= "$SNPs{$loc}{$indiv_1},$SNPs{$loc}{$indiv_2}";
print "\n\n***location of individual different SNPs (snps uniq to indiv)\n";
print "loc\tindividuals\tref\tsnps\n";
foreach my $loc (sort keys %diff_bwt_indiv){
  foreach my $strains (sort keys %{$diff_bwt_indiv{$loc}}){
    my $snps = $diff_bwt_indiv{$loc}{$strains}; 
    my ($chr,$pos) = $loc =~ /(\w+)\.(\d+)/;
    print "$loc\t$strains\t$ref{$loc}\t$snps\n";
  }
}

print "\n\n***location of individual different SNPs: homo only (snps uniq to indiv)\n";
print "loc\tindividuals\tref\tsnps\n";
my %count;
foreach my $loc (sort keys %diff_bwt_indiv){
  my $strain;
  foreach my $strains (sort keys %{$diff_bwt_indiv{$loc}}){
    my $snps = $diff_bwt_indiv{$loc}{$strains}; 
    #skip position if one snp is a het
    my ($chr,$pos) = $loc =~ /(\w+)\.(\d+)/;
    my $ref = $ref{"$chr.$pos"};
    next if $snps !~ /(.)\/\1,(.)\/\2/;   
    my @strains = split ',', $strains;
    my @snps = split ',' , $snps;
    for (my $i=0 ; $i<@strains ; $i++ ){
      my ($a1 , $a2) = $snps[$i] =~ /(.)\/(.)/;
      if ($ref ne $a1){
        $strain = $strains[$i];
      }
    }
   
    $count{$strain}++;
    print "$loc\t$strains\t$ref\t$snps\n";
  }
}
#print "\n\ncount of uniques Ref=T A/A,A/T is not counted for either strain:\n";
#foreach my $strain (keys %count){
#  print "$strain\t$count{$strain}\n";
#}

my %uniqHomo;
my %uniqHet;
my %mixed;
my $unq_hets;
my %total_hom;
my $total_shared_hom;
foreach my $loc (sort keys %binary){
  my @strains = keys %{$binary{$loc}};
  my @codes;
  ##0=same as ref
  ##1=homo diff from ref
  ##2=het 
  ##9=not enough data ./.
  foreach my $strain (@strains){
  #foreach my $strain (keys %{$binary{$loc}}){
    my $code = $binary{$loc}{$strain};
    push @codes ,$code;
  }
  my $codes = join '' , @codes;
  my $no_snp_count = $codes =~ tr/0/0/;
  my $homo_diff_from_ref_count = $codes =~ tr/1/1/;
  my $hets_count = $codes =~ tr/2/2/;
  my $missing_data_count = $codes =~ tr/9/9/;
 
  ## for total homo snps in each strain
  ## if all strains are a 1 or a zero, then if it is a 1 add up its total hom count;
  if ($hets_count == 0 and $missing_data_count == 0 and $homo_diff_from_ref_count > 0){
    for (my $i=0 ; $i < @codes ; $i++){
      my $code = $codes[$i];
      if ($code == 1){
        $total_hom{$strains[$i]}++;
      }
    }
  }
  ## for total shared homo snps
  ## if all strains are a 1, then if it is a 1 add up its total shared hom count;
  if ($hets_count == 0 and $missing_data_count == 0 and $homo_diff_from_ref_count == scalar @strains){
        $total_shared_hom++;
  }

  ## uniq homozygous snps
  ## if found only 1 homo diff from ref --> uniq homo snp
  if ($homo_diff_from_ref_count == 1 and $hets_count ==0 and $missing_data_count ==0){
   for (my $i=0 ; $i < @codes ; $i++){
     my $code = $codes[$i];
     if ($code == 1){
       $uniqHomo{$strains[$i]}++;
     }
   }
  }
  ## if found only 1 het diff from ref --> uniq het snp
  if ($homo_diff_from_ref_count == 0 and $hets_count ==1 and $missing_data_count ==0){
    for (my $i=0 ; $i < @codes ; $i++){
      my $code = $codes[$i];
       if ($code == 2){
         $uniqHet{$strains[$i]}++;
         my @snps_at_loc;
         my @strains = sort keys %{$SNPs{$loc}};
         foreach my $strain (@strains){
           push @snps_at_loc , $SNPs{$loc}{$strain};
         }
         my $snps_at_loc = join ',' , @snps_at_loc ;
         my $strains =  join ',' , @strains;
         $unq_hets .= "$loc\tref=$ref{$loc}\t$strains\t$snps_at_loc\n";
       }
    }
  }
  ## 2 strains both have SNPs, one HOMO, one HETE. could be Ref=T A/A,A/T. all other strains need to be 0
  if ($homo_diff_from_ref_count == 1 and $hets_count ==1 and $missing_data_count ==0){
    for (my $i=0 ; $i < @codes ; $i++){
      my $code = $codes[$i];
      $mixed{$loc}{$strains[$i]}=$code;
    }
  }
}
print "\n\ncount of unique HOMO SNPs. Ref=T A/A,T/T is counted. Ref=T A/A,A/T is not counted for either strain:\n";
foreach my $strain (keys %uniqHomo){
  my $count = $uniqHomo{$strain};
  print "$strain\t$count\n";
}
print "\n\ncount of total HOMO SNPs. Ref=T A/A,T/T is counted. Ref=T A/A,A/T is not counted for either strain:\n";
foreach my $strain (keys %total_hom){
  my $count = $total_hom{$strain};
  print "$strain\t$count\n";
}
print "\n\ncount of shared HOMO SNPs. Ref=T T/T,T/T is counted. Ref=T A/A,A/T is not counted for either strain:\n";
print "@strains\t$total_shared_hom\n";

print "\n\ncount of unique HET SNPs. Ref=T T/T,A/T is counted. Ref=T A/A,A/T is not counted for either strain:\n";
foreach my $strain (keys %uniqHet){
  my $count = $uniqHet{$strain};
  print "$strain\t$count\n";
}
print "**Locations unique HET SNPs.\n";
print $unq_hets ,"\n\n";


print "\n\ncount of locations with  2 strains with SNPs, one HOM one HET SNPs. Ref=T A/A,A/T is counted. Ref=T T/T,A/T is not counted for either strain:\n";
my $mixed_count=0;
my %mixed_strains;
foreach my $loc (keys %mixed){
  my @mixed_strains;
  $mixed_count++;
  foreach my $strain (keys %{$mixed{$loc}}){
    my $code = $mixed{$loc}{$strain};
    if ($code == 1 or $code == 2){
      push @mixed_strains, $strain;
    }
  }
  my @sorted = sort @mixed_strains;
  my $strains = join ',' , @sorted;
  $mixed_strains{$strains}++;
}
foreach my $strains(keys %mixed_strains){
  my $count = $mixed_strains{$strains};
  print "$strains\t$count\n";
}
print "$mixed_count\n";
