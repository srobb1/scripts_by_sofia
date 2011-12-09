#!/usr/bin/perl -w
use strict;
use File::Basename;
use Bio::Index::Fastq;


########Get User input############################################
my $blat_file  = shift;
my $fq_file_1  = shift;
my $len_cutoff = shift;
$len_cutoff = defined $len_cutoff ? $len_cutoff : 10;
##################################################################



########INDEX FASTQ###############################################
open INFQ, $fq_file_1 or die $!;
my ( $filename, $directories, $suffix ) = fileparse( $blat_file, qr/\.[^.]*/ );
my $fq_index = index_fastq($fq_file_1);
##################################################################



########OPEN OUTPUT FILES#########################################
open OUTMPING5, ">$filename.mping_five_prime.fa";
open OUTMPING3, ">$filename.mping_three_prime.fa";
open OUTFQ, ">$filename.mpingContainingReads$suffix";
##################################################################



####### PARSE BLAT ###############################################
####### STORE RESULTS ############################################
####### TRIM SEQS ################################################
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
    my $tStart   = $blat[15];
    my $tEnd     = $blat[16] - 1; #get all values into 1st base = 0 postion notation
    my $tLen     = $blat[14];
    if ( $strand eq '-' ) {
	my ($temp_tStart, $temp_tEnd) = ($tStart,$tEnd);
    	$tStart = $tLen - $temp_tEnd - 1;
    	$tEnd   = $tLen - $temp_tStart - 1;
    }
    my $addRecord;
    ##multi hits for same query read: mping match smaller therefore non-mping seq longer
    if ( exists $coord{$qName} ) {
        if ( $match < $coord{$qName}{match} ) {
            $addRecord = 1;
        }else {
            $addRecord = 0;
        }
    } else {
        $addRecord = 1;
    }
    if ($addRecord) {
	my $seq_obj = get_fq_record($qName, $fq_index);
	$coord{$qName}{seqObj}    = $seq_obj;
        $coord{$qName}{match}    = $match;
        $coord{$qName}{len}      = $qLen;
        $coord{$qName}{start}    = $qStart;
        $coord{$qName}{end}      = $qEnd;
        $coord{$qName}{tLen}     = $tLen;
        $coord{$qName}{mismatch} = $mismatch;
        $coord{$qName}{strand}   = $strand;
        if ( $strand eq '+' ) {
            $coord{$qName}{tStart} = $tStart;
            $coord{$qName}{tEnd}   = $tEnd;
        } else {
            $coord{$qName}{tStart} = $tLen - $tEnd - 1;
            $coord{$qName}{tEnd}   = $tLen - $tStart - 1;
        }
	trim_seqs(\%coord, $qName);
    }
}
##################################################################


sub trim_seqs{

    #want to cut and keep anything not matching to database TE
    my ( $coord, $qName ) = shift;

    my $seq_obj  = ${ $coord{$qName}{seqObj} };
    my $match    = ${ $coord{$qName}{match} };
    my $qLen     = ${ $coord{$qName}{len} };
    my $qStart   = ${ $coord{$qName}{start} };
    my $qEnd     = ${ $coord{$qName}{end} };
    my $tLen     = ${ $coord{$qName}{tLen} };
    my $mismatch = ${ $coord{$qName}{mismatch} };
    my $strand   = ${ $coord{$qName}{strand} };
    my $tStart   = ${ $coord{$qName}{tStart} };
    my $tEnd     = ${ $coord{$qName}{tEnd} };

    my ( $trimmed_seq, $trimmed_qual );
    my $header     = $seq_obj->id;
    my $seq        = $seq_obj->seq;
    my $qualHeader = '+';
    my $qual       = $seq_obj->qual_text;

    #query read overlaps 5' end of database TE & trimmed seq > cutoff
    if (    $tStart == 0
        and ( $qStart > $len_cutoff )
        and ( $mismatch == 0 )
      ) ##keep for later: allow a few mismatches## or ($mismatch / $len ) < 0.1 ) )
    {
        $trimmed_seq  = substr $seq,  0, $qStart;
        $trimmed_qual = substr $qual, 0, $qStart;
        my ( $tS, $tE, $qS, $qE ) =
          ( $tStart + 1, $tEnd + 1, $qStart + 1, $qEnd + 1 );
        print OUTMPING5
          ">$header $qS..$qE matches mping:$tS..$tE mismatches:$mismatch\n",
          substr( $seq, $qStart, ( $qEnd - $qStart + 1 ) ), "\n";
    }

    #query read overlaps 3' end of database TE & trimmed seq > cutoff
    elsif ( ( $tEnd == $tLen - 1 )
        and ( ( $qLen - $qEnd ) > $len_cutoff )
        and ( ( $mismatch == 0 ) )
      ) ##keep for later: allow a few mismatches## or ($mismatch / $len ) < 0.1 ) ) or ( $mismatch / $len < 0.1 ) ) )
    {
        $trimmed_seq  = substr $seq,  $qEnd + 1;
        $trimmed_qual = substr $qual, $qEnd + 1;
        my ( $tS, $tE, $qS, $qE ) =
          ( $tStart + 1, $tEnd + 1, $qStart + 1, $qEnd + 1 );
        print OUTMPING3
          ">$header $qS..$qE matches mping:$tS..$tE mismatches:$mismatch\n",
          substr( $seq, $qStart, ( $qEnd - $qStart + 1 ) ), "\n";
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

sub index_fastq {
    my $fq_file = shift;
    my ( $fname, $dir, $ext ) = fileparse( $fq_file, qr/\.[^.]*/ );
    my $index_file_name = $fname . ".index";
    my $inx             = Bio::Index::Fastq->new(
        '-filename'   => $index_file_name,
        '-write_flag' => 1
    );
    return $inx->make_index($fq_file);
}
sub get_fq_record{
	my $id = shift;
	my $index = shift;
	return $index->get_Seq_by_id($id);
	#returns Bio::Seq Obj
}
