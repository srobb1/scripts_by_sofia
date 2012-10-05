#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::Seq;

my $file  = shift;
my $db = shift; ##genome fasta

if ( !defined $file or !defined $db ) {
  die "./script file_with_range_info_to_be_collected fasta_file_of_genome";
}

#my $db_obj = Bio::DB::Fasta->new($fasta);
my %ref;

open IN, "$file" or die "Can't open $file $!";
my $fasta = $file;
$fasta =~ s/txt/fa/;
my @path = split /\// , $file;
my $name = pop @path;
my ($pre,$ext) = $name =~ /(\S+)\.(\S+)$/;
my $path = join '/' , @path;
open FAOUT, ">$path/flanking.$pre.fa";

my $count;
while ( my $line = <IN> ) {
  next if $line =~ /Not found/;
  chomp $line;
  #ZM_CACTA_93|methyl|chr10:1541752..1541998|chr10:1541752..1542063|complete|chr10:1541752..1542063
  my ($name,$methyl,$ref,$start,$end) = $line =~/^(\S+)\|(methyl.+)\|complete\|(\S+):(\d+)\.\.(\d+)$/;
  my $five_prime_start   = $start - 25;
  my $five_prime_end     = $start + 25 - 1;
  my $three_prime_start = $end - 25  + 1;
  my $three_prime_end   = $end + 25;
  my $five_prime      = getSeq_fastacmd( $ref, $five_prime_start,      $five_prime_end );
  my $three_prime   = getSeq_fastacmd( $ref, $three_prime_start,   $three_prime_end );
  print FAOUT ">$name|$methyl|complete|$ref:$start..$end|five_prime|$ref:$five_prime_start..$five_prime_end\n$five_prime\n";
  print FAOUT ">$name|$methyl|complete|$ref:$start..$end|three_prime|$ref:$three_prime_start..$three_prime_end\n$three_prime\n";

}

sub getSeq_fastacmd {
  my $ref   = shift;
  my $start = shift;
  my $end   = shift;
  my $seq;
  if (!exists $ref{$ref}){ 
    my $record = `fastacmd -d $db -s $ref`;
    my ($header, @seq) = split /\n/ , $record; 
    $seq = join '', @seq;
    $ref{$ref}=$seq;
  }else{
    $seq = $ref{$ref};
  
  }
  if ($seq){
    my $length = length($seq);
    $end = $length < $end ? $length : $end;
    if ( defined $seq ) {
      ##correct for 0 1 base notation
      $start--;
      $end--;
      my $subseq = substr ($seq,$start, ($end-$start+1));
      return $subseq;
    }
    else {
      return "subseq retrival issue";
    }
  }else {
    return "seqObj retrieval issue";
  }

}
