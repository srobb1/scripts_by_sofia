#!/usr/bin/perl -w
use strict;
use File::Basename;

my $blat_file  = defined $ARGV[0] ? $ARGV[0] : die "get_mping_fa_from_each_read.pl <blatfile_1> <fq_1>\n" ;
my $fq_file_1  = defined $ARGV[1] ? $ARGV[1] : die "get_mping_fa_from_each_read.pl <blatfile_1> <fq_1>\n" ;

open INBLAT, $blat_file, or die "Please provide a blat output file\n";

<INBLAT>;    #get rid of blat header lines
<INBLAT>;
<INBLAT>;
<INBLAT>;
<INBLAT>;

my %coord;
while ( my $line = <INBLAT> ) {
    chomp $line;
    my @blat = split /\t/, $line;

    my $mismatch = $blat[1];
    my $strand   = $blat[8];
    my $qName    = $blat[9];
    my $qLen     = $blat[10];
    my $qStart   = $blat[11];
    my $qEnd     = $blat[12];
    my $tStart   = $blat[15];
    my $tEnd     = $blat[16];
    my $tLen     = $blat[14];

    $coord{$qName}{len}      = $qLen;
    $coord{$qName}{start}    = $qStart;
    $coord{$qName}{end}      = $qEnd;
    $coord{$qName}{tLen}     = $tLen;
    $coord{$qName}{mismatch} = $mismatch;
    $coord{$qName}{strand}   = $strand;
    if ( $strand eq '+' ) {
        $coord{$qName}{tStart} = $tStart;
        $coord{$qName}{tEnd}   = $tEnd;
    }
    else {
        $coord{$qName}{tStart} = $tLen - $tEnd - 1;
        $coord{$qName}{tEnd}   = $tLen - $tStart - 1;
    }
}

open INFQ, $fq_file_1 or die $!;
my ( $filename, $directories, $suffix ) = fileparse( $fq_file_1, qr/\.[^.]*/ );
open OUTMPINGFA,     ">$filename.mping.from.reads.fa";
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
        my ( $tS, $tE, $qS, $qE ) = ( $tStart + 1, $tEnd + 1, $start + 1, $end + 1 );

        print OUTMPINGFA ">$header $qS..$qE matches mping:$tS..$tE mismatches:$mismatch\n",
              substr( $seq, $start, ( $end - $start + 1 ) ), "\n";
	}
}

