#!/usr/bin/perl -w
use strict;

my $blat_file          = shift;
my $db          = shift;

open (INBLAT, $blat_file) or die "Please provide a blat output file\n";

<INBLAT>;    #get rid of blat header lines
<INBLAT>;
<INBLAT>;
<INBLAT>;
<INBLAT>;

while ( my $line = <INBLAT> ) {
  chomp $line;
  my @blat     = split /\t/, $line;
  my $match    = $blat[0];
  my $mismatch = $blat[1];
  my $strand   = $blat[8];
  my $qName    = $blat[9];
  my $qLen     = $blat[10];
  my $qStart   = $blat[11];
  my $qEnd   = $blat[12] - 1; #get all values into 1st base = 0 postion notation
  my $tName   = $blat[13];
  my $tLen   = $blat[14];
  my $tStart = $blat[15];
  my $tEnd   = $blat[16] - 1; #get all values into 1st base = 0 postion notation
  my $block_qStarts = $blat[19];
  my ($block_qStart) = split ',', $block_qStarts;
  my $addRecord = 0;
  next if ($match+$mismatch) == $tLen;
  if ($tName =~ /pong_5prime/ and $match >= 28 and $tStart == 0 and (($qStart <= $qLen - 28))) { #and ($qStart >= 10) )){
#72      0       0       0       0       0       0       0       -       p00.FC89_1_p1.te_pong.ContainingReads1:1:2657:15304:Y:start     100     0       72      pong_5prime     120     0       72      1       72,     28,     0,
    $addRecord = 1;
  ##pong_3prime should overlap with end of seq
  }elsif ($tName =~ /pong_3prime/ and $match >= 18 and (($qStart <= $tLen - 18) and $tEnd == ($tLen-1))) {
    $addRecord = 1;
  }
  if ($addRecord) {
     #print "$tName\n$qName\n";
     my $seqRec = `fastacmd -d $db -s $qName`;
     die if !defined $seqRec;
     my ($header, @seq) = split /\n/ , $seqRec;
     my $seq = join '' , @seq;
     my $trim_seq;
     my $trimmed = 0;
     if ($qStart == 0){
       next if ($qEnd+1) >= $qLen;
       $trim_seq = substr $seq, $qEnd+1;
       $trimmed=1;
#       print "header:$header qN:$qName qS:$qStart qE:$qEnd qL:$qLen sL:" , length $seq ."\n" if !defined $trim_seq;
     }elsif($qEnd == ($qLen-1)){
       $trim_seq = substr $seq, 0 , $qStart;
       $trimmed=1;
     }else{
       $trimmed = 0;
       $trim_seq = "$qStart $qEnd $tLen" . ' notTrimmed';
     }
     if ($strand eq '-'){
       my $rev = reverse $seq;
       $seq = $rev;
       $seq  =~ tr /ATGCN/TACGN/;

       my $rev_trim = reverse $trim_seq;
       $trim_seq = $rev_trim;
       $trim_seq  =~ tr /ATGCN/TACGN/;
   
     }
     #print ">$tName.$qName\n$seq\n" if length $trim_seq > 10;
     print ">$tName.$qName.trim\n$trim_seq\n" if length $trim_seq > 10 and $trimmed;
  }
}


__END__
- strand
>1:1:2657:15304:Y:start
A123_0.te_containing_fq.fasta-TTGCTCTTTTGCTGACTTGTCACTGTATTAAATGCGCATGACACTCAAATGAAACACCCATTGTGACTGGCC TAAGAGTGAGAAACGAACAGTCAAGTAG
