#!/usr/bin/perl -w
use strict;
use Bio::SearchIO;

my $file = shift;
my $searchIO_obj = new Bio::SearchIO(-format => 'blast',  
                           -file   => $file);
my %hits;
while( my $result = $searchIO_obj->next_result ) {
  my $count = 0;
  my $q_name = $result->query_name;
  my $q_len = $result->query_length;
  $hits{$q_name}{len}=$q_len;
  while( my $hit = $result->next_hit ) {
    $count++;    
    my $evalue = $hit->significance;
    if ($evalue < 1e-10){
      my $desc = $hit->description;
      my $h_name = $hit->name;
      my $h_len = $hit->length;

      if ($desc =~ /gene/){
        $hits{$q_name}{hit}{$h_name}{desc}=$desc;
        $hits{$q_name}{hit}{$h_name}{evalue}=$evalue;
        $hits{$q_name}{hit}{$h_name}{h_len}=$h_len;
      }
    }
    last if $count == 5;
  }
}
print "q_name\tq_len\tevalue\thit_name\thit_len\tdesc\n";
foreach my $q (sort keys %hits){
  my $len = $hits{$q}{len};
  if (exists $hits{$q}{hit}){
    foreach my $h (sort { $hits{$q}{hit}{$a}{evalue} <=> $hits{$q}{hit}{$b}{evalue} } keys %{$hits{$q}{hit}}){
      my $evalue = $hits{$q}{hit}{$h}{evalue};
      my $desc = $hits{$q}{hit}{$h}{desc};
      my $len = $hits{$q}{hit}{$h}{h_len};
      print "$q\t$len\t$evalue\t$h\t$len\t$desc\n";
    }
  }else{
      print "$q\t$len\t-\t-\t-\n";
  }
  print "\n";
}
