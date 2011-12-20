#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;
use File::Basename;

## trims and filters reads, aligns to genome, splits by target.
## produces fq, sam and bam files for each individual target 

my $dir = '.';
my $genomeFasta;
my ( $minLength, $minQuality, $minPercent ) = ( 50, 20, 80 );
my $Q            = 33;
my $insertLength = 500;
my $mate_1_id    = "_1";
my $mate_2_id    = "_2";
my $split = 0;
my $filter_trim = 1;
my $tempDir = "/scratch";
#my $tempDir = "/dev/shm";
#my $start_host = `hostname`;
#$start_host =~ s/\s+//g;


GetOptions(
    'd|dir:s'          => \$dir,
    '1|mate_1_id:s'    => \$mate_1_id,
    '2|mate_2_id:s'    => \$mate_2_id,
    'g|genomeFasta:s'  => \$genomeFasta,
    'l|minLength:i'    => \$minLength,
    'q|minQuality:i'   => \$minQuality,
    'p|minPercent:i'   => \$minPercent,
    's|quality:i'      => \$Q,
    'i|insertLength:i' => \$insertLength, 
    'x|split:i'	       => \$split,
    'f|filter_trim:i'  => \$filter_trim,
    't|tempDir:s'      => \$tempDir,
    'h|help'           => \&getHelp,
);

if ( !defined $genomeFasta ) {
    print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
    &getHelp();
}

sub getHelp () {
    print "
usage:
./raw_paired_reads_2_split_by_target.pl [-d fq_file_directory] [-1 mate_pair_file_1_id][-2 mate_pair_file_2_id][-l minLength] [-q minQuality] [-p minPercent] [-s quality_offset] [-i insert_size] [-h] 

options:
-d STR		directory of original raw fq files (.fq not .fastq) [.]
-g STR		genome fasta file path [no default]
-l INT		min length for fastq_quality_trimmer [50]
-q INT		min quality score for fastq_quality_trimmer [20]
-p INT		min percent for fastq_quality_filter [80]
-s INT		quality score offset type Sanger(33) or Illumina(64) [33]
-i INT		insert library length [500]
-1 STR		file containing mate 1 id (ex reads_1.fq) [_1]
-2 STR		file containing mate 2 id (ex reads_2.fq) [_2]
-x INT	        split fq file into smaller files (1,000,000/file) yes=1 no=0 [0]	
-f INT	        run fastq_quality_filter and fastq_quality_trimmer yes=1 no=0 [1]	
-t STR		location to create temp directories [/dev/shm]
-h 		this message
";

    exit 1;
}

unless ( -e $genomeFasta ) {
    &getHelp();
    die "$genomeFasta does not exist. Check file name.\n";
}

my $genome_path = File::Spec->rel2abs($genomeFasta);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

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
my $dir_path = File::Spec->rel2abs($dir);
##split into smaller files
if ($split){
	mkdir "$current_dir/split_by_number_fq", 0777 unless -d "$current_dir/split_by_number_fq";    
	opendir( DIR, $dir_path ) || die "$!";
	foreach my $file ( readdir(DIR) ) {
	    my ( $volume, $directories, $filename ) = File::Spec->splitpath($file);
	    next unless ( $filename =~ /((\S+)($mate_1_id|$mate_2_id))\.(fastq|fq)$/ );
	    my ( $filename_base, $sampleName, $pairID, $suffix ) = ( $1, $2, $3, $4 );
	    `~/bin/fastq_split.pl -s 1000000 -o split_by_number_fq/ $dir_path/$file`;
	}
}
if ($split){
	$dir_path = "$current_dir/split_by_number_fq";
}
opendir( DIR, $dir_path ) || die "Can't Open $dir_path $!";
my $fq_ext;
foreach my $file ( readdir(DIR) ) {
    my ( $volume, $directories, $filename ) = File::Spec->splitpath($file);
    next unless ( $filename =~ /((\S+)($mate_1_id|$mate_2_id))\.(fastq|fq)$/ );
    my ( $filename_base, $sampleName, $pairID, $suffix ) = ( $1, $2, $3, $4 );
    push @${ $files{$sampleName} }, $filename_base;
    $fq_ext = $suffix;
}

my $desc;
my $ext;
if ($filter_trim){
	$desc = ".trimmed.filtered";      
	$ext = ".fq";                 
}else{
	$desc = "";  
	$ext = ".fq";
}


foreach my $sample ( sort keys %files ) {
    open OUTFILE, ">$current_dir/$sample.sh";
    print OUTFILE "#!/bin/bash\n\n";
    my (
        @trim_filter, @clean,               @aln,
        @sam,         @split_sam_by_target, @sam2fq,
        @sam2bam,     @mergeBam,            @catfq
    );

    print OUTFILE "tmp_dir=`mktemp --tmpdir=$tempDir -d`\n";
    print OUTFILE "cd \$tmp_dir\n";

    #foreach single file write the trim and filter and the aln commands
    foreach my $file ( sort @${ $files{$sample} } ) {
	if ($filter_trim){
		push @trim_filter, 
"fastq_quality_trimmer -Q$Q -l $minLength -t $minQuality -i $dir_path/$file.$fq_ext |fastq_quality_filter -Q$Q -q $minQuality -p $minPercent -v -o \$tmp_dir/$file.trimmed.filtered.fq";
        	push @aln,
"bwa aln -t 10 $genome_path \$tmp_dir/$file$desc.matched$ext > \$tmp_dir/$file$desc.matched.sai";
#"bwa aln -t 10 $genome_path \$tmp_dir/$file$desc$ext.matched > \$tmp_dir/$file$desc.matched.sai";
	}else{
		print OUTFILE "ln -s $dir_path/$file.$fq_ext \$tmp_dir/.\n";
        	push @aln,
"bwa aln -t 10 -q 10 $genome_path \$tmp_dir/$file$desc.matched$ext > \$tmp_dir/$file$desc.matched.sai";
#"bwa aln -t 10 -q 10 $genome_path \$tmp_dir/$file$desc$ext.matched > \$tmp_dir/$file$desc.matched.sai";
	}

    }

    #foreach potential paired sample write the following commands
    my ( $pair1, $pair2 );
    my $pairs = @${ $files{$sample} };
    if ( $pairs == 2 ) {
        ( $pair1, $pair2 ) = sort @${ $files{$sample} };
        push @clean,
"~/bin/clean_pairs.pl -1 \$tmp_dir/$pair1$desc$ext -2 \$tmp_dir/$pair2$desc$ext > \$tmp_dir/$sample.unPaired.fq";
        ##after cleaning the 2 paired files a file of unPaired reads is generated
        ##run bwa aln on this file
        ##and bwa samse
        push @clean, "if [ -s \$tmp_dir/$sample.unPaired.fq ] ; then bwa aln -t 10 $genome_path \$tmp_dir/$sample.unPaired.fq > \$tmp_dir/$sample.unPaired.sai ; fi";
        push @clean, "if [ -e \$tmp_dir/$sample.unPaired.sai ] ; then bwa samse  $genome_path \$tmp_dir/$sample.unPaired.sai \$tmp_dir/$sample.unPaired.fq   > \$tmp_dir/$sample.unPaired.sam ; fi";
        push @split_sam_by_target,"if [ -e \$tmp_dir/$sample.unPaired.sam ] ; then ~/bin/splitSam_byTarget.pl -s \$tmp_dir/$sample.unPaired.sam ; fi";
        push @sam,
"bwa sampe -a $insertLength $genome_path \$tmp_dir/$pair1$desc.matched.sai \$tmp_dir/$pair2$desc.matched.sai \$tmp_dir/$pair1$desc.matched$ext \$tmp_dir/$pair2$desc.matched$ext  > \$tmp_dir/$sample.sam";
        push @split_sam_by_target,
          "~/bin/splitSam_byTarget.pl -s \$tmp_dir/$sample.sam";
        push @sam2fq,
"for i in `ls \$tmp_dir/split_by_target` ; do ~/bin/sam2fq.pl \$tmp_dir/split_by_target/\$i ; done";
        push @sam2bam,
"for i in `ls \$tmp_dir/split_by_target` ; do ~/bin/sam2bam.pl \$tmp_dir/split_by_target/\$i $genome_path $sample; done";
#        push @sam2bam,
#"for i in `ls \$tmp_dir/split_by_target` ; do ~/bin/sam2bam.pl \$tmp_dir/split_by_target/\$i $genome_path $sample.unPaired; done";
    }
    elsif ( $pairs == 1 ) {
        ( $pair1, $pair2 ) = sort @${ $files{$sample} };
        push @sam,
"bwa samse $genome_path \$tmp_dir/$pair1$desc.sai \$tmp_dir/$pair1$desc.fq   > \$tmp_dir/$sample.sam";
        push @split_sam_by_target,
          "~/bin/splitSam_byTarget.pl -s \$tmp_dir/$sample.sam";
    }
    else {
        warn
"error: $sample has $pairs. This sample should have 2 pairs or just 1.\n";
    }
    	
    foreach my $trim_filter (@trim_filter) {
        print OUTFILE "$trim_filter\n\n";
    }
    foreach my $clean (@clean) {
        print OUTFILE "$clean\n\n";
    }
    foreach my $aln (@aln) {
        print OUTFILE "$aln\n\n";
    }
    foreach my $sam (@sam) {
        print OUTFILE "$sam\n\n";
    }
    foreach my $line (@split_sam_by_target) {
        print OUTFILE "$line\n\n";
    }
    foreach my $line (@sam2fq) {
        print OUTFILE "$line\n\n";
    }
    foreach my $line (@sam2bam) {
        print OUTFILE "$line\n\n";
    }
    print OUTFILE "mkdir -p $current_dir/sam_for_all_reads\n";
    print OUTFILE "cp \$tmp_dir/*sam $current_dir/sam_for_all_reads\n";
    if ($filter_trim){
    	print OUTFILE "mkdir -p $current_dir/fq_split_by_number_filtered\n";
    	print OUTFILE "cp \$tmp_dir/*.matched.fq \$tmp_dir/*unPaired.fq $current_dir/fq_split_by_number_filtered\n";
    }
    print OUTFILE "mkdir -p $current_dir/fq_split_by_chromosome\n";
    print OUTFILE "mkdir -p $current_dir/sam_split_by_chromosome\n";
    print OUTFILE "mkdir -p $current_dir/bam_split_by_chromosome\n";
    print OUTFILE "for i in `ls \$tmp_dir/split_by_target` ; do mkdir -p $current_dir/fq_split_by_chromosome/\$i ; done \n";
    print OUTFILE "for i in `ls \$tmp_dir/split_by_target` ; do mkdir -p $current_dir/sam_split_by_chromosome/\$i ; done \n";
    print OUTFILE "for i in `ls \$tmp_dir/split_by_target` ; do mkdir -p $current_dir/bam_split_by_chromosome/\$i ; done\n";
    print OUTFILE "for i in `ls \$tmp_dir/split_by_target` ; do ~/bin/move_files.pl \$tmp_dir/split_by_target/\$i $current_dir; done\n";
    print OUTFILE "cd $current_dir\n";
    #print OUTFILE "host_remote=`hostname`\n";
    #print OUTFILE "if [ $start_host = \$host_remote ]; then cd $current_dir; rm -rf \$tmp_dir; else ssh \$host_remote \"rm -rf \$tmp_dir\"; fi\n";
    print OUTFILE  "rm -rf \$tmp_dir";
}

