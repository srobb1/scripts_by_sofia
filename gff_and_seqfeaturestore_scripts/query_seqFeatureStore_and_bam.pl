#!/usr/bin/perl
use warnings;
use strict;
use Bio::DB::SeqFeature::Store;
use Bio::DB::Sam;

#file: getReads.pl

#######################################################################
# 1. This script gets all gene features from a seqfeature:store database
# 2. then gets all the snoRNAs (or some other feature) that are contained
# 	within the start and stops of the gene feature
# 3. then it uses the locations of the snoRNA feature to find RNAseq
# 	reads in a Bam file that overlap the start and stop.
#
# see this error '[bam_header_read] EOF marker is absent.' but it does not
# SEEM to effect the data. 	
#
# ./getReads.l > overlapping_reads.fa
#
#########################################################################


#############################################################################
## make a sam object (contains all data in bam file and all methods needed to
## interact with bam file
#############################################################################
my $sam = Bio::DB::Sam->new(
    -bam =>
      "/Library/WebServer/Documents/gbrowse2/databases/human/TEST.sorted.bam",
    -fasta => "/Library/WebServer/Documents/gbrowse2/databases/human/chr21.fa",
);

#############################################################################
## make a Bio::DB::SeqFeature::Store object (contains info about  HUMANSOFIA 
## database and all the methods needed to interact interact with the database
#############################################################################
my $human_db = Bio::DB::SeqFeature::Store->new(
    -adaptor => 'DBI::mysql',
    -dsn     => 'dbi:mysql:HUMANSOFIA',
    -user    => 'root',
);



# queries DB for all features of type 'gene'
my @features_type_gene = $human_db->get_features_by_type('gene');
foreach my $feature (@features_type_gene) {
    my $gene_name  = $feature->name;
    my $gene_start = $feature->start;
    my $gene_end   = $feature->end;
    my $ref        = $feature->ref;

    my $gene_info = "$gene_name($ref:$gene_start..$gene_end)";

    my $gene_type = $feature->type;

    if ( $gene_type !~ /pseudo/i and $gene_type !~ /RNA/ ) {

	# queries DB for all features that are within the boundaries of the 
	# gene feature
        my @features_by_location = $human_db->get_features_by_location(
            -seq_id => $ref,
            -start  => $gene_start,
            -end    => $gene_end
        );
        foreach my $feature (@features_by_location) {
            my $type_by_location = $feature->type;

            if ( $type_by_location =~ /gene:snoRNA/i ) {
                my $name  = $feature->name;
                my $start = $feature->start;
                my $end   = $feature->end;
                my $sub_feature_info = "$name($ref:$start..$end)";


		# subroutine for retrieving reads from the bam file
		# that overlap the start and stop provided
                my @read_info = getOverlappingReads( $ref, $start );
                foreach my $read_info (@read_info) {
                    print ">$gene_info $sub_feature_info $read_info";
                }
                
		@read_info = getOverlappingReads( $ref, $end );
                foreach my $read_info (@read_info) {
                    print ">$gene_info $sub_feature_info $read_info";
                }
            }

        }
    }
}

sub getOverlappingReads {
    my ( $ref, $pos ) = @_;
    my @read_info;
    my @alignments = $sam->get_features_by_location(
        -seq_id => $ref,
        -start  => $pos,
        -end    => $pos
    );
    foreach my $alignment (@alignments) {
        my $length  = $alignment->length;
        my $f_start = $alignment->start;
        my $f_end   = $alignment->end;
        my $f_seq   = $alignment->query->dna;

        push @read_info, "read:$f_start..$f_end len=$length\n$f_seq\n";
    }
    return @read_info;
}

