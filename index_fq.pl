#!/usr/bin/perl -w
use strict;
use File::Basename;
use Bio::Index::Fastq;


########Get User input############################################
my @fq_files  = @ARGV;
##################################################################



########INDEX FASTQ###############################################
my $fq_index = index_fastq(@fq_files);
##################################################################


sub index_fastq {
    my $fq_files = @_;
    my ( $fname, $dir, $ext ) = fileparse( $fq_files[0], qr/\.[^.]*/ );
    my $index_file_name = $fname . ".index";
    my $inx             = Bio::Index::Fastq->new(
        '-filename'   => $index_file_name,
        '-write_flag' => 1
    );
    return $inx->make_index(@fq_files);
}
