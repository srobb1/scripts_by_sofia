#!/usr/bin/perl -w
use strict;
use Data::Dumper;

my $dir = shift;
##make options to change these
my $mate_file_1 = '_1.';
my $mate_file_2 = '_2.';
my $mate_file_unpaired = '_unPaired.';

my @files_1 = <$dir*$mate_file_1*fq*>;
my @files_2 = <$dir*$mate_file_2*fq*>;
my @files_unpaired = <$dir*$mate_file_unpaired*fq*>;


my %flanking_fq;

for (my $i = 0 ; $i < @files_1 ; $i++){
  my $file_1 = $files_1[$i];
  $file_1 =~ s/$mate_file_1//;
  for (my $j = 0 ; $j < @files_2 ; $j++){
    my $file_2 = $files_2[$j];
    $file_2 =~ s/$mate_file_2//;
    if ($file_1 eq $file_2){
      $flanking_fq{$file_1}{1}=$files_1[$i];
      $flanking_fq{$file_1}{2}=$files_2[$j];
      if (@files_unpaired){ 
        for (my $k = 0 ; $k < @files_unpaired ; $k++){
          my $file_unpaired = $files_unpaired[$k];
          if ($file_1 eq $file_unpaired){
            $flanking_fq{$file_1}{unpaired}=$files_unpaired[$k];
            last;
          }
        }
      }
      last; #if $file_1 eq $file_2 and we are finished with unpaired go back to $i loop
    }
  }
}
print Dumper \%flanking_fq;
