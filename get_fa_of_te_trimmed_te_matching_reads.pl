#!/usr/bin/perl -w
use strict;
use File::Basename;
use Bio::SeqIO;
my $blat_file  = shift;
my $fa_file_1  = shift;
my $len_cutoff = shift;
my $mismatch_allowance = shift; 
$len_cutoff = defined $len_cutoff ? $len_cutoff : 10;
$mismatch_allowance = defined $mismatch_allowance ? $mismatch_allowance : 0.1;

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

#open INFA, $fa_file_1 or die $!;
my ( $filename, $directories, $suffix ) = fileparse( $blat_file, qr/\.[^.]*/ );
my $TE = "unspecified";
my $FA = "unspecified";
#$fa.te_$TE.blatout
if ($filename =~ /(\S+)\.te_(\S+)/){
	$FA = $1;
	$TE = $2;
}
my $seqIO_obj = Bio::SeqIO->new(-file => "$fa_file_1", -format => 'fasta');
while ( my $seqObj = $seqIO_obj->next_seq ) {
    my $header = $seqObj->id;
    my $seq        = $seqObj->seq;

    if ( exists $coord{$header} ) {
        my $start    = $coord{$header}{start};
        my $len      = $coord{$header}{len};
        my $end      = $coord{$header}{end};
        my $tStart   = $coord{$header}{tStart};
        my $tEnd     = $coord{$header}{tEnd};
        my $tLen     = $coord{$header}{tLen};
        my $mismatch = $coord{$header}{mismatch};
        my $match = $coord{$header}{match};
	my $strand   = $coord{$header}{strand};
        #want to cut and keep anything not matching to database TE
        my ( $trimmed_seq, $trimmed_qual );

        #query read overlaps 5' end of database TE & trimmed seq > cutoff
        if (    $tStart == 0
            and ( ($len-($match+$mismatch)) > $len_cutoff )
            and ($mismatch / $len ) <= $mismatch_allowance  
          ) 
        {
	    my ( $tS, $tE, $qS, $qE ) =
              ( $tStart + 1, $tEnd + 1, $start + 1, $end + 1 );
	    if ($strand eq '-'){
		## start at the end of the match and go to end of string
                $trimmed_seq  = substr $seq,  $end+1; 
	    }else{ ## strand is positive 
                $trimmed_seq  = substr $seq,  0, $start;
	    }
	    next if length $trimmed_seq < $len_cutoff; 
        }

        #query read overlaps 3' end of database TE & trimmed seq > cutoff
        elsif ( ( $tEnd == $tLen - 1 )
            and  ($len-($match+$mismatch) > $len_cutoff )
            and ($mismatch / $len ) <= $mismatch_allowance  
          ) ##keep for later: allow a few mismatches## or ($mismatch / $len ) < 0.1 ) ) or ( $mismatch / $len < 0.1 ) ) )
        {
            my ( $tS, $tE, $qS, $qE ) =
              ( $tStart + 1, $tEnd + 1, $start + 1, $end + 1 );
	      if ($strand eq '-'){
		    $trimmed_seq  = substr $seq,  0 , $start;
	      }else { ## strand is +
		    $trimmed_seq  = substr $seq,  $end + 1;
	      }
	    next if length $trimmed_seq < $len_cutoff; 
        }
        if ( defined $trimmed_seq ) {
            print ">$header\n$trimmed_seq\n";
        }
    }
}

