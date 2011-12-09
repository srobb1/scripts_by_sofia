#!/usr/bin/perl -w

use strict;
use File::Spec;

my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);
my $genome_path = '/home_stajichlab/robb/wesslerlab-shared/Rice/Genome/index/osa1r6_plus_Mito_Chloroplast.fa';
#my $genome_path = '/home_stajichlab/robb/wesslerlab-shared/Rice/Genome/index/IRGSP_B5_plusMito_plusPlastid.fa;
my $dir = shift;
my $dir_path = File::Spec->rel2abs($dir);
#$dir_path .= '/';
opendir( DIR, $dir ) || die "$!";
foreach my $file ( readdir(DIR) ) {
    #my ($volume,$directories,$filename) = File::Spec->splitpath( $file );
    my $sample;
    if ($file =~ /(\S+)\.sam$/){
	$sample = $1;
	open OUTSH, ">$dir_path/$sample.scratch.sam2bam.sh" or die "$dir_path/$sample.scratch.sam2bam.sh ".$!;
	print OUTSH "#!/bin/sh\n\n";
	print OUTSH "tmp_dir=\`mktemp --tmpdir=/scratch -d\`\n";
	print OUTSH 'cd $tmp_dir',"\n";	
	print OUTSH "samtools view -b -S -T $genome_path  $dir_path/$sample.sam -o \$tmp_dir/$sample.bam\n" ;
	print OUTSH "samtools sort  \$tmp_dir/$sample.bam  \$tmp_dir/$sample.sorted\n";
	print OUTSH "samtools index  \$tmp_dir/$sample.sorted.bam\n";
	print OUTSH "cd $dir_path\n";
	print OUTSH 'remote_host=`hostname`' , "\n";
	print OUTSH "echo \"scp \$remote_host:\$tmp_dir/* $dir_path/.\" > $dir_path/$sample.toMove.sh\n";
	print OUTSH "echo \"ssh \$remote_host rm -rf \$tmp_dir\" >> $dir_path/$sample.toMove.sh\n";
	close OUTSH;
     }

}
