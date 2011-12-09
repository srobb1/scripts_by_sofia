#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;

## creates a shell script to cat fq reads into files 
## based on chromosome 

my $dir = '.';
my $pair_id_a = "_1";
my $pair_id_b = "_2";
my $unpaired_id = "unPaired";
my $fq_ext = "fq";
my $file_out_base;
GetOptions(
    'd|dir:s'          => \$dir,
    'a|pair_id_a:s'    => \$pair_id_a,
    'b|pair_id_b:s'    => \$pair_id_b,
    'u|unpaired_id:s'  => \$unpaired_id,		
    'x|fq_ext:s'       => \$fq_ext,
    'o|file_out_base:s'=> \$file_out_base,
    'h|help'           => \&getHelp,
	
);

if (! defined $file_out_base){

	print "Please supply a base name for output files\n";
	print "Example:	strain.chromosome01 would result in
			strain.chromosome01_1.fq
			strain.chromosome01_2.fq\n";
	&getHelp;
}

sub getHelp () {
    print "
usage:
./SCRIPTNAME_HERE.pl [-d fq_file_directory] 

options:
-d STR		directory of fq files (.fq not .fastq) [.]
-a STR		identifier of the left or first read pairs [_1]
-b STR		identifier of the right or second read pairs [_2]
-u STR		identifier of the unpaired reads
-x STR		file extension of fastq files [fq]
-h 		this message
";

    exit 1;
}

my $dir_path = File::Spec->rel2abs($dir);
opendir( DIR, $dir ) || die "$!";


my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);
print "
--------------------
Parameters Used:
current directory $current_dir
-d  $dir_path
-a  $pair_id_a
-b  $pair_id_b
-u  $unpaired_id
-o  $file_out_base
-x  $fq_ext

Run the resulting shell script to retrieve your concatenated files.

\n
---------------------
\n";

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $uniqID = $sec.$min.$hour;
my @unpaired;
my @paired_1;
my @paired_2;
my %files;
foreach my $file ( readdir(DIR) ) {
    my ( $volume, $directories, $filename ) = File::Spec->splitpath($file);

    next unless ( $filename =~ /((\S+?)($pair_id_a|$pair_id_b))\.$fq_ext$/ );
    my ( $filename_base, $sampleName, $pairID ) = ( $1, $2, $3 );
  
    my $complete_file_path = "$dir_path/$filename";


    if ($filename =~ /$unpaired_id/){
	push @unpaired, $complete_file_path;
    }elsif ($filename =~ /$pair_id_a.$fq_ext/){
	push @paired_1, $complete_file_path;
    }elsif ($filename =~ /$pair_id_b.$fq_ext/){
  	push @paired_2, $complete_file_path;
    }	

	
}

my $pair_1 = join " ", @paired_1;
my $pair_2 = join " ", @paired_2;
my $unpaired = join " ", @unpaired;


open OUTFILE, ">$current_dir/$uniqID.$file_out_base.cat_fq.sh";
print OUTFILE "#!/bin/bash\n\n";
print OUTFILE 'tmp_dir=`mktemp --tmpdir=/scratch -d`', "\n";
print OUTFILE "cd \$tmp_dir\n";
print OUTFILE "cat $pair_1 > \$tmp_dir/$file_out_base"."$pair_id_a".".$fq_ext\n";
print OUTFILE "cat $pair_2 > \$tmp_dir/$file_out_base"."$pair_id_b".".$fq_ext\n";
print OUTFILE "cat $unpaired > \$tmp_dir/$file_out_base".".$unpaired_id".".$fq_ext\n";
print OUTFILE "cd $current_dir\n";
print OUTFILE 'remote_host=`hostname`' , "\n";
print OUTFILE "echo \"scp \$remote_host:\$tmp_dir/* $dir_path/.\" >> $current_dir/$uniqID.$file_out_base.toMove.sh\n";
print OUTFILE "echo \"ssh \$remote_host rm -rf \$tmp_dir\" >> $current_dir/$uniqID.$file_out_base.toMove.sh\n";
print OUTFILE "echo \"Done!!\"\n";

