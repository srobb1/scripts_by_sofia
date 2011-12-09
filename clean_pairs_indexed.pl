#!/usr/bin/perl -w

use strict;
use File::Spec;
use Getopt::Long;
use lib '/home_stajichlab/robb/bin/';
use Bio::Index::Fastq;

##takes a sam file and splits it into separate files based on the
##targets sequences.

##$file_1 and _2 should be mate fq files
my $file_1;
my $file_2;



GetOptions(
    '1|file_1:s' => \$file_1,
    '2|file_2:s' => \$file_2,
    'h|help'  => \&getHelp,
);


sub getHelp () {
    print "
This script takes 2 paired fastq files and matches mates and prints unPaired reads to STDOUT. 
Please name the fastq files with .fq extention.
usage:
./cleanUp_Pairs.pl [-1 fastq file 1] [-2 fastq file 2][-h] 

options:
-1 STR          sam file 1 [required, no default]
-2 STR          sam file 2 [required, no default]
-h              this message
";
    exit 1;
}

if ( !defined $file_1 and !defined $file_2 ) {
    print "\n\nMust provide a two files containing paired reads\n\n";
    &getHelp;
    exit 1;
}

my $file_1_path = File::Spec->rel2abs($file_1);
my $file_2_path = File::Spec->rel2abs($file_2);


my $fq_1_inx = Bio::Index::Fastq->new('-filename' => "index.1", '-write_flag' => 1);
my $fq_2_inx = Bio::Index::Fastq->new('-filename' => "index.2", '-write_flag' => 1);
$fq_1_inx->make_index($file_1_path);
$fq_2_inx->make_index($file_2_path);
my $sample;
##sample_1.fq
##sampe_pair1.fq
if ($file_1 =~ /(\S+?)_?(\S+)\.(fq|fastq)/){
	$sample = $1;
}

my %pairs;

my $out_seqIO_1 = Bio::SeqIO->new('-format' => 'Fastq','-filename' => "$file_1.matched");
my $out_seqIO_2 = Bio::SeqIO->new('-format' => 'Fastq','-filename' => "$file_2.matched");
my $out_seqIO_unMatched = Bio::SeqIO->new('-format' => 'Fastq','-fh' => \*STDOUT);
my @ids = $fq_1_inx->get_all_primary_ids;
print "len = scalar @ids\n";
foreach my $header ( @ids){
	print "**$header\n";
   my $header_to_store;
   my $seqObj = $fq_1_inx->get_Seq_by_id($header);
   if ($header =~ /(\S+)(\.r|\.f)/){
	$header_to_store = $1;
    }elsif($header =~ /(\S+)(\/1|\/2)/){
	$header_to_store = $1;
    }else {
	$header_to_store = $header;
    }
    $pairs{$header_to_store} = $seqObj;
	
}
foreach my $header ( $fq_2_inx->get_all_ids){
    my $header_to_store;
    if ($header =~ /(\S+)(\.r|\.f)/){
        $header_to_store = $1;
    }elsif($header =~ /(\S+)(\/1|\/2)/){
        $header_to_store = $1;
    }else {
        $header_to_store = $header;
    }    
    my $seqObj_2 = $fq_2_inx->get_Seq_by_id($header_to_store);

    my $seqObj_1 = $pairs{$header_to_store};
    if ( exists $pairs{$header_to_store} ) {
        $out_seqIO_1->write_seq($seqObj_1);
        $out_seqIO_2->write_seq($seqObj_2);

	#remove the found keys, so that only the unparied 
	#from the first file are in hash at end of loop
        delete $pairs{$header_to_store};
    }
    else {    #print any seq in file2 that was not in file one
    	$out_seqIO_unMatched->write_seq($seqObj_2);
    }
}

#will print any remaing and unpair seqs from the first file into NOMATCH
foreach my $header ( sort keys %pairs ) {
    my $seqObj = $pairs{$header};
    $out_seqIO_unMatched->write_seq($seqObj);
}
