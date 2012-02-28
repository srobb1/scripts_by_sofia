#!/usr/bin/perl -w
use strict;
use File::Spec;
#combines all the fq files in one directory
#provide the direcotry name
## for multiple directories try this:
## for i in `ls` ; do cat_fq_shell.pl $i prefix; done
if (!defined @ARGV){
  die "run command like this: for i in `ls` ; do cat_fq_shell.pl \$i prefix [/tmp_dir default:/sractch] [clean 1|0 default:1]; done\n";
}

my $dir = shift;
my $prefix = shift;
my $tempDir = shift;
my $clean = shift;
$prefix = !defined $prefix ? '' : $prefix . '.';
$clean = !defined $clean ? 1 : 0;
my $dir_path = File::Spec->rel2abs($dir);

$tempDir = defined $tempDir ? $tempDir : '/scratch';
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

my @dirs = split '/' , $dir_path;
my $lowest_dir = pop @dirs;
my $one_up = join '/' , @dirs;

open SH, ">$current_dir/$lowest_dir.cat_fq.sh";
print SH "#!/bin/bash\n\n";

print SH "echo mktemp\n";
print SH "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
print SH "echo \"tmp_dir=\$tmp_dir\"\n";
print SH "cd \$tmp_dir\n";

my $mate_1 = "$lowest_dir"."_1.fq" ; 
my $mate_2 = "$lowest_dir"."_2.fq" ; 
my $unpaired = "$lowest_dir"."_unPaired.fq"; 

print SH "cat $dir_path/*_1.fq > \$tmp_dir/$mate_1\n"; 
print SH "cat $dir_path/*_2.fq > \$tmp_dir/$mate_2\n"; 
print SH "if [ -s $dir_path/$unpaired ] ; then cat $dir_path/*unPaired.fq >  \$tmp_dir/$unpaired.tmp ; fi\n";
#print SH "cat $dir_path/*unPaired.fq >  \$tmp_dir/$unpaired.tmp\n";

if ($clean){
  print SH "clean_pairs.pl -1 \$tmp_dir/$mate_1 -2 \$tmp_dir/$mate_2 > \$tmp_dir/$unpaired.tmp2\n";
  print SH "if [ -e \$tmp_dir/$unpaired.tmp ] ; then cat \$tmp_dir/$unpaired.tmp \$tmp_dir/$unpaired.tmp2 > \$tmp_dir/$unpaired ; else mv \$tmp_dir/$unpaired.tmp2 \$tmp_dir/$unpaired ; fi\n";
  $mate_1 = "$lowest_dir"."_1.matched.fq"; 
  $mate_2 = "$lowest_dir"."_2.matched.fq"; 
}
print SH "cp \$tmp_dir/$mate_1 $one_up/$prefix$mate_1\n";
print SH "cp \$tmp_dir/$mate_2 $one_up/$prefix$mate_2\n";
print SH "if [ -s \$tmp_dir/$unpaired  ] ; then  cp \$tmp_dir/$unpaired $one_up/$prefix$unpaired ; fi\n";

print SH "cd $current_dir\n";
print SH "rm -rf \$tmp_dir\n";
