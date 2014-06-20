#!/usr/bin/perl -w

use strict;
use File::Spec;
#for i in `seq 1 12` ; do ~/bin/sam2bam_scratch.pl genome_path Chr$i ;done
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);
my $genome_path = shift; 
my $dir = shift;
my $dir_path = File::Spec->rel2abs($dir);
opendir( DIR, $dir ) || die "$!";
foreach my $file ( readdir(DIR) ) {
    my $sample;
    if ($file =~ /(\S+)\.sam$/){
	$sample = $1;
	open OUTSH, ">$current_dir/$sample.scratch.sam2bam.sh" or die "$current_dir/$sample.scratch.sam2bam.sh ".$!;
	print OUTSH "#!/bin/sh\n\n";
	print OUTSH "tmp_dir=\`mktemp --tmpdir=/scratch -d\`\n";
	print OUTSH 'cd $tmp_dir',"\n";	
	print OUTSH "samtools view -b -S -T $genome_path  $dir_path/$sample.sam -o \$tmp_dir/$sample.bam\n" ;
	print OUTSH "samtools sort  \$tmp_dir/$sample.bam  \$tmp_dir/$sample.sorted\n";
	print OUTSH "samtools index  \$tmp_dir/$sample.sorted.bam\n";
	print OUTSH "cd $dir_path\n";
	print OUTSH "cp \$tmp_dir/* $dir_path/.\n";
	print OUTSH "rm -rf \$tmp_dir\n";
	close OUTSH;
     }

}
