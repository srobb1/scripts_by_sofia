#!/usr/bin/perl -w

use strict;
use File::Spec;

###java -jar /home/robb/stajichlab-shared/pkg/Picard/picard/SamToFastq.jar INPUT=p00.FC52_7.sam FASTQ=test.pair1.fq SECOND_END_FASTQ=test.pair2.fq VALIDATION_STRINGENCY=LENIENT

my $dir = shift;
my $dir_path = File::Spec->rel2abs($dir);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

opendir( DIR, $dir ) || die "$!";
foreach my $file ( readdir(DIR) ) {
    my $filename = "$dir_path/$file";
    if ($file =~ /^(\S+)\.sam$/){
    	my $filebase = $1;
    	open OUTSH, ">$dir_path/$filebase.sam2fq.sh"; 
	print OUTSH "#!/bin/bash\n\n";
	print OUTSH "java -jar /home_stajichlab/robb/stajichlab-shared/pkg/Picard/picard/SamToFastq.jar INPUT=$dir_path/$file FASTQ=$dir_path/$filebase.fq VALIDATION_STRINGENCY=LENIENT\n";
	print OUTSH 'grep -A3 -P \.f$  '."$dir_path/$filebase.fq | grep -v".' ^--$ >'." $dir_path/$filebase"."_f.fq\n"; 
	print OUTSH 'grep -A3 -P \.r$  '."$dir_path/$filebase.fq | grep -v".' ^--$ >'." $dir_path/$filebase"."_r.fq\n"; 
	close OUTSH;    
    }
}
