#!/usr/bin/perl -w
use strict;
use File::Spec;
#combines all the fq files in one directory
#provide the direcotry name
## for multiple directories try this:
## for i in `ls` ; do cat_fq_shell.pl $i prefix; done
my $dir = shift;
my $prefix = shift;
$prefix = !defined $prefix ? '' : $prefix . '.';
my $dir_path = File::Spec->rel2abs($dir);

my $tempDir = defined $ARGV[2] ? $ARGV[2]  : '/scratch';
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

my @dirs = split '/' , $dir_path;
my $lowest_dir = pop @dirs;
my $one_up = join '/' , @dirs;

open SH, ">$current_dir/$lowest_dir.cat_fq.sh";
print SH "#!/bin/bash\n\n";

print SH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
print SH "cd \$tmp_dir\n";

my $mate_1 = "$lowest_dir"."_1.fq" ; 
my $mate_2 = "$lowest_dir"."_2.fq" ; 
my $unpaired = "$lowest_dir"."_unPaired.fq"; 

print SH "cat $dir_path/*_1.fq > \$tmp_dir/$mate_1\n"; 
print SH "cat $dir_path/*_2.fq > \$tmp_dir/$mate_2\n"; 
print SH "cat $dir_path/*unPaired.fq >  \$tmp_dir/$unpaired\n";

print SH "cp \$tmp_dir/$mate_1 $one_up/$prefix$mate_1\n";
print SH "cp \$tmp_dir/$mate_2 $one_up/$prefix$mate_2\n";
print SH "if [ -s \$tmp_dir/$unpaired  ] ; then  cp \$tmp_dir/$unpaired $one_up/$prefix$unpaired ; fi\n";

print SH "cd $current_dir\n";
print SH "rm -rf \$tmp_dir\n";
