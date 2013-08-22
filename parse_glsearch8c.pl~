#!/usr/bin/perl -w
use strict;

my $file = shift; ##glsearch c8 output

#19.2.1  H.3.5   100.00  61      0       0       1       61      1       61      8.6e-83 48.9

open IN , $file;

my %hits;
my %matches;

while (my $line = <IN> ){
  next if $line =~ /^#/;
  chomp $line ;
  my ($query, $hit, $id) = split /\s+/ , $line;
  next if $id < 80;
  $hits{$query}{$hit} = $id;
  next if $query eq $hit ;  
## should this be a '.' or a '\t'??
  if (exists $matches{"$query\t$hit"}){
    my $stored_id = $matches{"$query\t$hit"};
    if ($id > $stored_id){
      $matches{"$query\t$hit"} = $id;
    }
  }elsif (exists $matches{"$hit\t$query"}){
    my $stored_id = $matches{"$hit\t$query"};
    if ($id > $stored_id){
      $matches{"$hit\t$query"} = $id;
    }
  }
  else{
    $matches{"$query\t$hit"} = $id;
  }
}


my @header = sort keys %hits ;
print "\t" , join "\t" , @header , "\n";
foreach my $query (@header){
  my $row = "$query\t";
  foreach my $hit (@header){
    my $id = exists $hits{$query}{$hit} ?  $hits{$query}{$hit} : '-' ;
    $row .= "$id\t";
  }
  $row =~ s/\t$/\n/;
  print $row;
}

print "\n\nmatch1\tmatch2\t%id\n";
foreach my $match (sort {$matches{$b} <=> $matches{$a}}keys %matches){
  print "$match\t$matches{$match}\n";
}
