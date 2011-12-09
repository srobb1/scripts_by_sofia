#!/usr/bin/perl -w
##modified 11-3-2011
use strict;
use File::Spec;
#provide a directory of bam files to merge
# ex: a directory for one chromosome
# for i in `seq 1 12` ; do ~/bin/merge_bam_and_create_merged_sam_shell.pl  Chr$i ; done
#  or
# for i in `ls` ; do ~/bin/merge_bam_and_create_merged_sam_shell.pl  $i ; done
my $dir = shift;
my $prefix = shift;
$prefix = !defined $prefix ? '' : $prefix.'.'; 
my $dir_path = File::Spec->rel2abs($dir);
my $tempDir = defined $ARGV[2] ? $ARGV[2]  : '/scratch';
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

my $start_host = `hostname`;
$start_host =~ s/\s+//g;

my @dirs = split '/' , $dir_path;
my $lowest_dir = pop @dirs;
my $one_up = join '/' , @dirs;
pop @dirs;
my $two_up = join '/', @dirs;

open SH, ">$current_dir/$lowest_dir.mergeBam.sh";
print SH "#!/bin/bash\n\n";

print SH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
print SH "cd \$tmp_dir\n";



#merge bam files
print SH "samtools merge \$tmp_dir/$lowest_dir.merged.bam $dir_path/*bam\n";
#`samtools merge -f $one_up/$lowest_dir.merged.bam $dir_path/*bam`;

#sort and index merged bam
print SH "samtools sort  \$tmp_dir/$lowest_dir.merged.bam  \$tmp_dir/$lowest_dir.merged.sorted\n";
#`samtools sort  $one_up/$lowest_dir.merged.bam  $one_up/$lowest_dir.merged.sorted`;

print SH "samtools index  \$tmp_dir/$lowest_dir.merged.sorted.bam\n";
#`samtools index  $one_up/$lowest_dir.merged.sorted.bam`;


#convert merged and sorted bam to a sam
print SH "samtools view -h \$tmp_dir/$lowest_dir.merged.sorted.bam -o \$tmp_dir/$lowest_dir.merged.sorted.sam\n";
#`samtools view -h $one_up/$lowest_dir.merged.sorted.bam -o $two_up/sam_split_by_chromosome/$lowest_dir.merged.sorted.sam`;

print SH "for i in `ls *.sam*` ; do cp \$tmp_dir/$i $two_up/sam_split_by_chromosome/$prefix$i ; done\n";
print SH "for i in `ls *.bam*` ; do cp \$tmp_dir/$i $one_up/$prefix$i ; done\n";
print SH "cd $current_dir\n";
print SH "rm -rf \$tmp_dir\n";
##while the script is running it is on the remote node so i do not have to ssh to rm the files
#print SH "host_remote=`hostname`\n";
#print SH "if [ $start_host = \$host_remote ]; then cd $current_dir; rm -rf \$tmp_dir; else ssh \$host_remote \"rm -rf \$tmp_dir\"; fi\n";
