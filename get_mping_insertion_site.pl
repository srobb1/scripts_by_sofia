#!/usr/bin/perl -w
use strict;

my $sorted_bam  = shift;
my $usr_target  = shift;
my $genome_path = shift;
my $genome_seq;
my $id;
##get chromosome sequence for substr of flanking seq
open GENOME, "$genome_path" or die "Cannot open genome fasta: $genome_path $!";
my $i = 0;
while ( my $line = <GENOME> ) {
    last if $i > 1;
    chomp $line;
    if ( $line =~ />(\S+)/ ) {
        $i++;
        $id = $1;
        if ( $id ne $usr_target ) {
            warn "-t $usr_target is not found in genome fasta.";
            warn
              "$id is found in genome fasta and will be used to parse results.";
            $usr_target = $id;
        }
    }
    else {
        $genome_seq .= $line;
    }
}

my @sorted_sam = `samtools view $sorted_bam`;

my $last_start = 0;
my $last_end   = 0;
my %mpingInsertions;
my $count = 0;
my @bin   = (0);
foreach my $line (@sorted_sam) {
    chomp $line;
    my (
        $name,  $flag, $target, $start, $four,
        $cigar, $six,  $seven,  $eight, $seq
    ) = split /\t/, $line;
    my $len = length $seq;
    my $end = $len + $start - 1;
    next if $target ne $usr_target;

    ## test to see if we are still within one insertion event or a different one
    ## is this seq aligned to same region, is it in range
    my $padded_start = $bin[0] - 5;
    my $padded_end   = $bin[-1] + 5;
    if (   ( $start >= $padded_start and $start <= $padded_end )
        or ( $end >= $padded_start and $end <= $padded_end ) )
    {
        push @bin, $start, $end;
        @bin = sort @bin;
        my $first_3bp = uc( substr $seq, 0, 3 );
        my $last_3bp = uc( substr $seq, -3 );
        TAA_TTA_check( $count, $first_3bp, $start,   'left',  $name );
        TAA_TTA_check( $count, $last_3bp,  $end - 2, 'right', $name );
    }
    else {
        ## if start and end do not fall within last start and end
        ## we now have a different insertion event
        $count++;
        my $first_3bp = uc( substr $seq, 0, 3 );
        my $last_3bp = uc( substr $seq, -3 );
        TAA_TTA_check( $count, $first_3bp, $start,   'left',  $name );
        TAA_TTA_check( $count, $last_3bp,  $end - 2, 'right', $name );

        #reset last_start, last_end, @bin
        @bin        = ( $start, $end );
        $last_start = $start;
        $last_end   = $end;
    }
}

my $event = 0;
open OUTFASTA, ">$usr_target.mping_insertion_sites.fa"         or die $!;
open OUTALL,   ">$usr_target.mping_insertions.all.txt"         or die $!;
open OUTGFF,   ">$usr_target.mping_insertion_sites.gff"        or die $!;
open OUTTABLE, ">$usr_target.mping_insertion_sites.table.txt"  or die $!;
open OUTLIST,  ">$usr_target.mping_insertion_sites.reads.list" or die $!;
print OUTGFF "##gff-version	3\n";
##output for students in tab delimited table
print OUTTABLE
"chromosome\tinsertion_site\tleft_flanking_read_count\tright_flanking_read_count\tleft_flanking_seq\tright_flanking_seq\n";

foreach my $insertionEvent ( sort { $a <=> $b } keys %mpingInsertions ) {
    foreach my $start (
        sort { $a <=> $b }
        keys %{ $mpingInsertions{$insertionEvent} }
      )
    {
        my $start_count = $mpingInsertions{$insertionEvent}{$start}{count};
        my $left_count  = $mpingInsertions{$insertionEvent}{$start}{left};
        my $right_count = $mpingInsertions{$insertionEvent}{$start}{right};
        my @reads       = @{ $mpingInsertions{$insertionEvent}{$start}{reads} };

        if ( defined $left_count and defined $right_count ) {
            $event++;
            my $coor                  = $start + 2;
            my $flank_len             = 50;
            my $zero_base_coor        = $coor - 1;
            my $sub_string_start      = $zero_base_coor - $flank_len + 1;
            my $seq_start             = $coor - $flank_len + 1;
            my $seq_end               = $coor + $flank_len;
            my $left_flanking_ref_seq = substr $genome_seq, $sub_string_start,
              $flank_len;
            my $right_flanking_ref_seq = substr $genome_seq,
              $zero_base_coor + 1, $flank_len;
            print OUTTABLE
"$usr_target\t$coor\tleft_flanking_read_count=$left_count\tright_flanking_read_count=$right_count\tleft_flanking_seq=$left_flanking_ref_seq\tright_flanking_seq=$right_flanking_ref_seq\n";
            print OUTGFF
"$usr_target\t.\tinsertion_site\t$coor\t$coor\t.\t.\t.\tID=mping_insertion_site.$usr_target.$coor; left_flanking_read_count=$left_count; right_flanking_read_count=$right_count; left_flanking_seq=$left_flanking_ref_seq; right_flanking_seq=$right_flanking_ref_seq\n";
            print OUTFASTA
">$usr_target.$coor $usr_target:$seq_start..$seq_end\n$left_flanking_ref_seq$right_flanking_ref_seq\n";
            print OUTALL
"$usr_target\t$coor\tC:$start_count\tR:$right_count\tL:$left_count\n";
            print OUTLIST "$usr_target:$coor\t", join( ",", @reads ), "\n";
        }
        else {
            my $coor = $start + 2;
            $left_count  = defined $left_count  ? $left_count  : 0;
            $right_count = defined $right_count ? $right_count : 0;
            print OUTALL
"$usr_target\t$coor\tC:$start_count\tR:$right_count\tL:$left_count\n";
        }
    }
}
print OUTALL
  "C=total read count, R=right hand read count, L=left hand read count\n";
print OUTALL
"total confident insertions identified by a right AND left mping flanking sequence (C>1,R>0,L>0)= $event\n";

sub TAA_TTA_check {
    my ( $event, $seq, $start, $pos, $read ) = @_;
    if ( $seq eq 'TAA' or $seq eq 'TTA' ) {
        $mpingInsertions{$event}{$start}{count}++;
        $mpingInsertions{$event}{$start}{left}++  if $pos eq 'left';
        $mpingInsertions{$event}{$start}{right}++ if $pos eq 'right';
        push @{ $mpingInsertions{$event}{$start}{reads} }, $read;
    }
}
