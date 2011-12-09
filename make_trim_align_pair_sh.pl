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
my ( $minLength, $minQuality, $minPercent ) = ( 50, 20, 80 );
my $Q            = 33;
my $insertLength = 500;

GetOptions(
    'd|dir:s'          => \$dir,
    'g|genomeFasta:s'  => \$genomeFasta,
    'l|minLength:i'    => \$minLength,
    'q|minQuality:i'   => \$minQuality,
    'p|minPercent:i'   => \$minPercent,
    's|quality:i'      => \$Q,
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
./run_trim_align_pair.pl [-d fq_file_directory] [-l minLength] [-q minQuality] [-p minPercent] [-s quality_offset] [-i insert_size] [-h] 

options:
-d STR		directory of fq files (.fq not .fastq) [.]
-g STR		genome fasta file path [no default]
-l INT		min length for fastq_quality_trimmer [50]
-q INT		min quality score for fastq_quality_trimmer [20]
-p INT		min percent for fastq_quality_filter [80]
-s INT		quality score offset type Sanger(33) or Illumina(64) [33]
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
-l  $minLength
-q  $minQuality
-p  $minPercent
-s  $Q
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
   
    next unless ( $filename =~ /((\S+?)_(\S*pair)?(1|2))\.fq$/ );
    my ( $filename_base, $sampleName, $pairID ) = ( $1, $2, $4 );

    push @${ $files{$sampleName} }, $filename_base;

}

foreach my $sample ( sort keys %files ) {
    open OUTFILE, ">$current_dir/$sample.sh";
    print OUTFILE "#!/bin/bash\n\n";
    my ( @trim_filter, @clean, @aln, @sam );

    #foreach single file write the trim and filter and the aln commands
    foreach my $file ( sort @${ $files{$sample} } ) {
	push @trim_filter,
          "fastq_quality_trimmer -Q$Q -l $minLength -t $minQuality -i $dir_path"
          . "$file"
          . ".fq |fastq_quality_filter -Q$Q -q $minQuality -p $minPercent -v -o $current_dir/$file.trimmed.filtered.fq";
        push @aln,
"bwa aln -t 10 $genome_path $current_dir/$file.trimmed.filtered.fq.matched > $current_dir/$file.trimmed.filtered.matched.sai";

    }

    #foreach potential paired sample write the following commands
    my ( $pair1, $pair2 );
    my $pairs = @${ $files{$sample} };
    if ( $pairs == 2 ) {
        ( $pair1, $pair2 ) = sort @${ $files{$sample} };
        push @clean,
"~/bin/cleanUpPairs_fq.pl $current_dir/$pair1.trimmed.filtered.fq $current_dir/$pair2.trimmed.filtered.fq > $current_dir/$sample.unPaired.fq";
	##after cleaning a file of unPaired reads is generated
	##run bwa aln on this file
	##and bwa samse
        push @clean, "bwa aln -t 10 $genome_path $current_dir/$sample.unPaired.fq > $current_dir/$sample.unPaired.sai";
	push @clean,
"bwa samse  $genome_path $current_dir/$sample.unPaired.sai $current_dir/$sample.unPaired.fq   > $current_dir/$sample.unPaired.sam";
	push @sam,
"bwa sampe -a $insertLength $genome_path $current_dir/$pair1.trimmed.filtered.matched.sai $current_dir/$pair2.trimmed.filtered.matched.sai $current_dir/$pair1.trimmed.filtered.fq.matched $current_dir/$pair2.trimmed.filtered.fq.matched  > $current_dir/$sample.sam";
    }
    elsif ( $pairs == 1 ) {
        ( $pair1, $pair2 ) = sort @${ $files{$sample} };
        push @sam,
"bwa samse $genome_path $current_dir/$pair1.trimmed.filtered.sai $current_dir/$pair1.trimmed.filtered.fq   > $current_dir/$sample.sam";
    }
    else {
        warn
"error: $sample has $pairs. This sample should have 2 pairs or just 1.\n";
    }

    foreach my $trim_filter (  @trim_filter ) {
        print OUTFILE "echo \"$trim_filter\"\n";
        print OUTFILE "$trim_filter\n\n";
    }
    foreach my $clean ( @clean ) {
        print OUTFILE "echo \"$clean\"\n";
        print OUTFILE "$clean\n\n";
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

