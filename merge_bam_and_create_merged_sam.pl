#!/usr/bin/perl -w
use strict;
use File::Spec;
#provide a directory of bam files to merge
# ex: a directory for one chromosome
# for i in `seq 1 12` ; do ~/bin/merge_bam_and_create_merged_sam.pl  Chr$i ; done
my $dir = shift;
my $dir_path = File::Spec->rel2abs($dir);

my @dirs = split '/' , $dir_path;
my $lowest_dir = pop @dirs;
my $one_up = join '/' , @dirs;
pop @dirs;
my $two_up = join '/', @dirs;


#merge bam files
#print "samtools merge $one_up/$lowest_dir.merged.bam $dir_path/*bam\n";
`samtools merge -f $one_up/$lowest_dir.merged.bam $dir_path/*bam`;

#sort and index merged bam
#print "samtools sort  $one_up/$lowest_dir.merged.bam  $one_up/$lowest_dir.merged.sorted\n";
`samtools sort  $one_up/$lowest_dir.merged.bam  $one_up/$lowest_dir.merged.sorted`;
#print "samtools index  $one_up/$lowest_dir.merged.sorted.bam\n";
`samtools index  $one_up/$lowest_dir.merged.sorted.bam`;


#convert merged and sorted bam to a sam
#print "samtools view -h $one_up/$lowest_dir.merged.sorted.bam -o $two_up/sam_split_by_chromosome/$lowest_dir.merged.sorted.sam\n";
`samtools view -h $one_up/$lowest_dir.merged.sorted.bam -o $two_up/sam_split_by_chromosome/$lowest_dir.merged.sorted.sam`;

