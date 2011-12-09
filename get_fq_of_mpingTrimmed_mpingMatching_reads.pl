#!/usr/bin/perl -w
use strict;
use File::Basename;

my $blat_file  = shift;
my $fq_file_1  = shift;
my $len_cutoff = shift;
$len_cutoff = defined $len_cutoff ? $len_cutoff : 10;

open INBLAT, $blat_file, or die "Please provide a blat output file\n";

<INBLAT>;    #get rid of blat header lines
<INBLAT>;
<INBLAT>;
<INBLAT>;
<INBLAT>;

my %coord;
while ( my $line = <INBLAT> ) {
    chomp $line;
    my @blat     = split /\t/, $line;
    my $match    = $blat[0];
    my $mismatch = $blat[1];
    my $strand   = $blat[8];
    my $qName    = $blat[9];
    my $qLen     = $blat[10];
    my $qStart   = $blat[11];
    my $qEnd     = $blat[12] - 1; #get all values into 1st base = 0 postion notation
    my $tLen     = $blat[14];
    my $tStart   = $blat[15];
    my $tEnd     = $blat[16] - 1; #get all values into 1st base = 0 postion notation
    my $block_qStarts = $blat[19];
    my ($block_qStart) = split ',', $block_qStarts; 
    my $addRecord;
    ##multi hits for same query read: mping match smaller therefore non-mping seq longer
    if ( exists $coord{$qName} ) {
        if ( $match < $coord{$qName}{match} ) {
            $addRecord = 1;
        }
        else {
            $addRecord = 0;
        }
    }
    else {
        $addRecord = 1;
    }
    if ($addRecord) {
        $coord{$qName}{match}    = $match;
        $coord{$qName}{len}      = $qLen;
        $coord{$qName}{start}    = $qStart;
        $coord{$qName}{end}      = $qEnd;
        $coord{$qName}{tLen}     = $tLen;
        $coord{$qName}{mismatch} = $mismatch;
        $coord{$qName}{strand}   = $strand;
        $coord{$qName}{tStart} = $tStart;
        $coord{$qName}{tEnd}   = $tEnd;
        ##blat reports the coordinates in the positive strand for
	##both the query and the target. the coords of the revseq 
	##can be found in column19-qStarts for - strand query
	#if ( $strand eq '-' ) {
	#    $coord{$qName}{qStart} = $block_qStart;
	#    $coord{$qName}{qEnd} = $qLen - $qStart;
        #}
    }
}

open INFQ, $fq_file_1 or die $!;
my ( $filename, $directories, $suffix ) = fileparse( $blat_file, qr/\.[^.]*/ );
open OUTFQ,     ">$filename.mpingContainingReads$suffix";
open OUTMPING5, ">$filename.mping_five_prime.fa";
open OUTMPING3, ">$filename.mping_three_prime.fa";
while ( my $line = <INFQ> ) {
    chomp $line;
    my $header = $line;
    $header = substr( $header, 1 );
    chomp( my $seq        = <INFQ> );
    chomp( my $qualHeader = <INFQ> );
    chomp( my $qual       = <INFQ> );

    if ( exists $coord{$header} ) {
        my $start    = $coord{$header}{start};
        my $len      = $coord{$header}{len};
        my $end      = $coord{$header}{end};
        my $tStart   = $coord{$header}{tStart};
        my $tEnd     = $coord{$header}{tEnd};
        my $tLen     = $coord{$header}{tLen};
        my $mismatch = $coord{$header}{mismatch};
	my $strand   = $coord{$header}{strand};
        #want to cut and keep anything not matching to database TE
        my ( $trimmed_seq, $trimmed_qual );

        #query read overlaps 5' end of database TE & trimmed seq > cutoff
        if (    $tStart == 0
            and ( $start > $len_cutoff )
            and (  $mismatch == 0  )
          ) ##keep for later: allow a few mismatches## or ($mismatch / $len ) < 0.1 ) )
        {
            $trimmed_seq  = substr $seq,  0, $start;
            $trimmed_qual = substr $qual, 0, $start;
            my ( $tS, $tE, $qS, $qE ) =
              ( $tStart + 1, $tEnd + 1, $start + 1, $end + 1 );
              my $subseq = substr( $seq, $start, ( $end - $start + 1 ) );
	      if ($strand eq '-'){
		($subseq = reverse $subseq) =~ tr/AaGgTtCcNn/TtCcAaGgNn/;
	      }
	      print OUTMPING5
              ">$header $qS..$qE matches mping:$tS..$tE mismatches:$mismatch\n$subseq\n";
              #substr( $seq, $start, ( $end - $start + 1 ) ), "\n";
        }

        #query read overlaps 3' end of database TE & trimmed seq > cutoff
        elsif ( ( $tEnd == $tLen - 1 )
            and ( ( $len - $end ) > $len_cutoff )
            and ( ( $mismatch == 0 ) )
          ) ##keep for later: allow a few mismatches## or ($mismatch / $len ) < 0.1 ) ) or ( $mismatch / $len < 0.1 ) ) )
        {
            $trimmed_seq  = substr $seq,  $end + 1;
            $trimmed_qual = substr $qual, $end + 1;
            my ( $tS, $tE, $qS, $qE ) =
              ( $tStart + 1, $tEnd + 1, $start + 1, $end + 1 );
              my $subseq = substr( $seq, $start, ( $end - $start + 1 ) );
	      if ($strand eq '-'){
		($subseq = reverse $subseq) =~ tr/AaGgTtCcNn/TtCcAaGgNn/;
	      }
            print OUTMPING3
              ">$header $qS..$qE matches mping:$tS..$tE mismatches:$mismatch\n$subseq\n";
              #substr( $seq, $start, ( $end - $start + 1 ) ), "\n";
        }
        if ( defined $trimmed_seq ) {
            print "\@$header\n";
            print "$trimmed_seq\n";
            print "$qualHeader\n";
            print "$trimmed_qual\n";
        }
        ##mping containing reads
        ##any read that was in the blat file is written here
        print OUTFQ "\@$header\n";
        print OUTFQ "$seq\n";
        print OUTFQ "$qualHeader\n";
        print OUTFQ "$qual\n";
    }
}

