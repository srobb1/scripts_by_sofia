#!/usr/bin/perl -w
use strict;
#use inserts.characterixed.txt
#chromosome.pos  avg_flankers    spanners        status
#HEG4_2_nf       mping   TTA     Chr1.588834     7.5     29      new_insertion
#HEG4_2_nf       mping   TTA     Chr1.1132977    45.5    0       homozygous


my %inserts;
my @insert_files = @ARGV;
foreach my $file (@insert_files){
  open IN, "$file" or die "Can't open $file $!\n";
  while (my $line = <IN>){
   chomp $line;
   next if $line =~ /chromosome.pos.+avg_flan/;
      my ($strain, $te, $tsd, $pos, $spanners, $flankers, $class) = split /\t/ , $line;
      $class =~ s/\?//;
      $inserts{$pos}{$strain}=$class;
  }
}
my %uniq;
foreach my $pos (keys %inserts){
  my @strains = keys %{$inserts{$pos}};
  my $strains = join ( ',' , (sort @strains));
  $uniq{$strains}{count}++;
  my @classes;
  foreach my $strain (sort @strains){
    ## sort classes in same order as strains
    push @classes , $inserts{$pos}{$strain};
  }
  my $classes = join ( ',' , (@classes));
  $uniq{$strains}{pos}{$pos}=$classes;
}
my %class_count;
foreach my $strains (keys %uniq){
  foreach my $pos ( keys %{$uniq{$strains}{pos}}){
    my $class_count;
    my $classes  = $uniq{$strains}{pos}{$pos};
    $class_count{$strains}{$classes}++;;
    print "$strains\t$pos\t$classes\n";
  }
}
print "\n\n";
foreach my $strains (keys %class_count){
  foreach my $classes (keys %{$class_count{$strains}}){
    my $count = $class_count{$strains}{$classes};
    print "$strains\t$classes\t$count\n";
  }
}
print "\n\n";
foreach my $strains (keys %uniq){
  my $count = $uniq{$strains}{count};
  print "$strains\t$count\n";
}

