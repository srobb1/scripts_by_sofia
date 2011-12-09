#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;

## creates a shell script for each pair of illumina reads
## creates commands for:
## fastq_quality_trim
## fastq_quality_filter
## bwa align
## bwa sampe or samse

my $dir = '.';
my $genomeFasta;
my $insertLength = 500;

GetOptions(
    'd|dir:s'          => \$dir,
    'g|genomeFasta:s'  => \$genomeFasta,
    'i|insertLength:i' => \$insertLength,
    'h|help'           => \&getHelp,
);

if ( !defined $genomeFasta ) {
    &getHelp();
    die "\n\nPlease provide reference genome by using -g Genome fasta path\n";
}

sub getHelp () {
    print "
usage:
./make_align_sam.pl [-d fq_file_directory] [-i insert_size] [-h] 

options:
-d STR		directory of fq files (.fq not .fastq) [.]
-g STR		genome fasta file path [no default]
-i INT		insert library length [500]
-h 		this message
";

    exit 1;
}

my $dir_path = File::Spec->rel2abs($dir);
$dir_path .= '/';
opendir( DIR, $dir ) || die "$!";

unless ( -e $genomeFasta ) {
    &getHelp();
    die "$genomeFasta does not exist. Check file name.\n";
}

my $genome_path = File::Spec->rel2abs($genomeFasta);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);
print "
--------------------
Parameters Used:
current directory $current_dir
-d  $dir_path
-g $genome_path
-i  $insertLength
\n
---------------------
\n";

my $bwa_index = "bwa index -a bwtsw $genome_path";
unless ( -e "$genome_path.rsa" ) {
    open GINDEX, ">$current_dir/genome_indexing.sh";
    print GINDEX "#!/bin/bash\n\n";
    print GINDEX "echo \"$bwa_index\"\n";
    print GINDEX "$bwa_index\n";
    print GINDEX "chmod ago+rwx $genome_path*\n";
    close GINDEX;
}

my %files;
foreach my $file ( readdir(DIR) ) {
    my ($volume,$directories,$filename) = File::Spec->splitpath( $file );

   next unless ($filename =~ /(fq|fastq)(\.uniq)?$/);
    my ( $filename_base, $sampleName, $pairID );  
    if ($filename =~ /(\S+unPaired)\.fq\.uniq/){
	( $filename_base, $sampleName, $pairID ) = ( $1, $1, 1 )
    }
    elsif ( $filename =~ /((\S+?)_(\S*pair)?(1|2))\.fq\.uniq/ ){
    	( $filename_base, $sampleName, $pairID ) = ( $1, $2, $4 );
    }
    push @${ $files{$sampleName} }, $filename_base;

}

foreach my $sample ( sort keys %files ) {
    open OUTFILE, ">$current_dir/$sample.bwa-aln-sampe-samse.sh";
    print OUTFILE "#!/bin/bash\n\n";
    my (  @aln, @sam );

    #foreach single file write the trim and filter and the aln commands
    foreach my $file ( sort @${ $files{$sample} } ) {
        push @aln,
"bwa aln -t 10 $genome_path $current_dir/$file.fq.uniq > $current_dir/$file.uniq.sai";

    }

    #foreach potential paired sample write the following commands
    my ( $pair1, $pair2 );
    my $pairs = @${ $files{$sample} };
    if ( $pairs == 2 ) {
        ( $pair1, $pair2 ) = sort @${ $files{$sample} };
	push @sam,
"bwa sampe -a $insertLength $genome_path $current_dir/$pair1.uniq.sai $current_dir/$pair2.uniq.sai $current_dir/$pair1.fq.uniq $current_dir/$pair2.fq.uniq  > $current_dir/$sample.uniq.sam";
    }
    elsif ( $pairs == 1 ) {
        ( $pair1, $pair2 ) = sort @${ $files{$sample} };
        push @sam,
"bwa samse $genome_path $current_dir/$pair1.uniq.sai $current_dir/$pair1.fq.uniq   > $current_dir/$sample.uniq.sam";
    }
    else {
        warn
"error: $sample has $pairs. This sample should have 2 pairs or just 1.\n";
    }

    foreach my $aln ( @aln ) {
        print OUTFILE "echo \"$aln\"\n";
        print OUTFILE "$aln\n\n";
    }
    foreach my $sam (  @sam ) {
        print OUTFILE "echo \"$sam\"\n";
        print OUTFILE "$sam\n\n";
    }
    print OUTFILE "echo \"Done!!\"\n";
}

