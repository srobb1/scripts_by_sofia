#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use Cwd;
use strict;
if ( !defined @ARGV ) {
    &getHelp();
}

my $genomeFasta = 'NONE';
my $te_fasta;
my $target             = 'NONE';
my $len_cutoff         = 10;
my $mismatch_allowance = 0;
my $fq_dir;
my $exper              = 'not.given';
my $mate_file_1        = '_1.';
my $mate_file_2        = '_2.';
my $mate_file_unpaired = '_unPaired.';
my $outdir;

GetOptions(
    'e|exper:s'       => \$exper,
    'o|outdir:s'      => \$outdir,
    'd|fq_dir:s'      => \$fq_dir,
    'g|genomeFasta:s' => \$genomeFasta,
    't|te_fasta:s'    => \$te_fasta,
    'l|len_cutoff:i'  => \$len_cutoff,
    'm|mismatch:f'    => \$mismatch_allowance,
    '1|mate_1_id:s'   => \$mate_file_1,
    '2|mate_2_id:s'   => \$mate_file_2,
    'u|unpaired_id:s' => \$mate_file_unpaired,
    'h|help'          => \&getHelp,
);
my $current_dir;

if ( defined $outdir and -d $outdir ) {
    $current_dir = File::Spec->rel2abs($outdir);
}
else {
    $current_dir = cwd();
}
my $mapping = 1;
if ( !defined $genomeFasta ) {
    print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
    &getHelp();
}
elsif ( $genomeFasta eq 'NONE' ) {
    print
"You did not provide a genome fasta, if you proceed only reads containing the TE will be found, no mapping of insertions will be performed\n";
    print "Proceed without mapping?\n";
    my $answer;
    while ( $answer = <STDIN> ) {

        # Exit if it was just spaces (or just an enter)
        last if $answer =~ /^\s*|\n$/;
    }
    if ( $answer =~ /n/i ) {
        &getHelp();
    }
    else {
        $mapping = 0;
    }
}
elsif ( !-e $genomeFasta ) {
    print "$genomeFasta does not exist. Check file name.\n";
    &getHelp();
}
else {
    my $count = 0;
    my $seq;
    open( INFASTA, "$genomeFasta" ) || die "$!\n";
    while ( my $line = <INFASTA> ) {
        if ( $line =~ /^>(\S+)/ ) {
            $count++;
            last if $count > 1;
            $target = $1;
            $seq .= $line;
        }
        else {
            $seq .= $line;
        }
    }
    close INFASTA;

    if ( $count > 1 ) {
        open NEWGENOMEFASTA, "$current_dir/$target.fa" or die $!;
        print NEWGENOMEFASTA $seq;
        close NEWGENOMEFASTA;
        print
"\n\nGenome Fasta file:$genomeFasta contains more than one sequence.\n";
        print
"\n\nOnly the first sequence, $target, will be used to map TE insertions\n\n";
        $genomeFasta = "$current_dir/$target.fa";
    }
}
if ( !defined $te_fasta ) {
    print
"\n\nPlease provide fasta file containing transposable elements by using -t TE fasta path\n";
    &getHelp();
}
elsif ( !-e $te_fasta ) {
    print "$te_fasta does not exist. Check file name.\n";
    &getHelp();
}
else {
    my $first_line = `head -n1 $te_fasta`;
    if ( $first_line !~ /^>\S+\s+\S+/ ) {
        die
"The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME TSD\nSEQUENCE\n";
    }
}
my @fq_files;
my %fq_files;
if ( !defined $fq_dir ) {
    print "\n\nPlease provide a directory of paired fastq files\n";
    &getHelp();
}
elsif ( !-d $fq_dir ) {
    print
"\n\nCheck the spelling or location of $fq_dir, Please provide a directory of paired fastq files\n";
    &getHelp();
}
else {
    my $fq_path = File::Spec->rel2abs($fq_dir);
    @fq_files = <$fq_path/*fq>;
    my @fastq_files = <$fq_path/*fastq>;
    push @fq_files, @fastq_files;
    if ( scalar @fq_files == 0 ) {
        print "Must provide at least 1 short read file\n";
        &getHelp();
    }
}

sub getHelp {
    print "
usage:
./find_TE_insertionSites.pl [-t TE_fasta_file][-g chromosome_genome_fasta][-d dir_of_fq][-e short_sample_name][-h] 

options:

**required:
-g STR          single chromosome genome fasta file path [no default]
-t STR          fasta containing 1 or more nucleotide sequences of transposable elements with TSD in the desc [no default]
-d STR          directory of paired and unpaired fastq files (paired _1.fq & _2.fq) (.fq or .fastq is acceptable)  [no default]

**recommended:
-e STR          Short Sample name, will be used in the output files to create IDs for the insert (ex. A123) [not.given]

**optional:
-o STR          output directory, needs to exist, will not create, full path [cwd] 
-l INT          len cutoff for the te trimmed reads to be aligned [10] 
-m FRACTION     mismatch allowance for alignment to TE (int, ex 0.1) [0] 
-1 STR		string to identify mate 1 paired files [_1.]
-2 STR          string to identify mate 2 paired files [_2.]
-u STR          string to identify unpaied files [_unPaired.] 
-h              this message

SAMPLE TE FASTA
>mping TTA
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATG
ATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTT
TCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGT
CCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAA
CTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGT
TTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATT
GTGACTGGCC

";

    exit 1;
}

my $genome_path = File::Spec->rel2abs($genomeFasta) if $mapping;
my $te_path = File::Spec->rel2abs($te_fasta);

my @genome_dirs = split '/', $genome_path;
my $genome_file = pop @genome_dirs;
my $genome_dir  = join '/', @genome_dirs;
if ( !-e "$genome_dir/$genome_file.bowtie_build_index.1.ebwt" and $mapping ) {
    `bowtie-build -f $genome_path $genome_dir/$genome_file.bowtie_build_index`;
}

#create an index of genome fasta
`samtools faidx $genome_path`;

my @fq;
my @fa;

#foreach my $ref_fq_files ( keys %fq_files ) {

foreach my $fq (@fq_files) {

    #foreach my $fq ( @{ $fq_files{$ref_fq_files} } ) {
    my $fq_path = File::Spec->rel2abs($fq);
    push @fq, $fq_path;
    my $fa = $fq;
    if ( $fa =~ s/\.(fq|fastq)$/.fa/ ) {
        push @fa, $fa;
        if ( !-e $fa ) {
            open INFQ,  $fq_path or die $!;
            open OUTFA, ">$fa"   or die $!;

            #print "converting fq to fa $fq -> $fa\n";
            while ( my $header = <INFQ> ) {
                my $seq         = <INFQ>;
                my $qual_header = <INFQ>;
                my $qual        = <INFQ>;

                die "ERROR: expected \'\@\' but saw $header"
                  if substr( $header, 0, 1 ) ne '@';

                print OUTFA ">", substr( $header, 1 );
                print OUTFA $seq;
            }
            close INFQ;
            close OUTFA;
        }
    }
    else {
        print
"$fq does not seem to be a fastq based on the file extension. It should be fq or fastq\n";
        &getHelp();
    }
}

#}

#split TE fasta into single record fastas
my $num = 0;
open( INFASTA, "$te_fasta" ) || die "$!\n";
while ( my $line = <INFASTA> ) {
    if ( $line =~ /^>/ ) {
        $num++;
    }
}
close(INFASTA);
my @te_fastas;
my %TSD;

#put in a new directory structure workingDir/date-te-search/tefilename/all-newly-created-files
#create new te fasta file
my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst ) =
  localtime(time);
$mday = ( length $mday == 1 ) ? "0$mday" : $mday;
$mon++;
$mon = ( length $mon == 1 ) ? "0$mon" : $mon;
$year += 1900;

my $top_dir = $mday . $mon . $year . "_teSearch";
open( INFASTA, "$te_fasta" ) || die "$!\n";
my $i = 0;
while ( my $line = <INFASTA> ) {
    if ( $line =~ /^>(\S+)\s+(\S+)/ ) {
        my $id = $1;
        $TSD{$id} = $2;
        if ( $i > 0 ) {
            close(OUTFASTA);
            $i = 0;
        }
        my $te_file = "$id.fa";
        $te_file =~ s/\|/_/g;
        my $te_dir = "$current_dir/$top_dir/$id";
        push @te_fastas, "$te_dir/$te_file";
        `mkdir -p $te_dir`;
        open( OUTFASTA, ">$te_dir/$te_file" ) or die "$!\n";
        print OUTFASTA $line;
        $i++;
    }
    else {
        print OUTFASTA $line;
    }
}
close(INFASTA);
close(OUTFASTA);

#foreach TE fasta blat against target chromosome and parse and find insertion sites
foreach my $te_path (@te_fastas) {
    my @path     = split '/', $te_path;
    my $te_fasta = pop @path;
    my $path     = join '/', @path;
    my $TE       = $te_fasta;
    $TE =~ s/\.fa//;

    #blat fa files against te.fa
    my @flanking_fq;
    my @flanking_fq_mates;
    my $fq_file_count = scalar @fq;
    for ( my $i = 0 ; $i < $fq_file_count ; $i++ ) {
        my $fa = $fa[$i];
        my $fq = $fq[$i];

        #remove and save filename part of path
        my @fa_path = split '/', $fa;
        my $fa_name = pop @fa_path;
        $fa_name =~ s/\.fa$//;

`blat -minScore=10 -tileSize=7 $path/$te_fasta $fa $path/$fa_name.te_$TE.blatout`
          if !-e "$path/$fa_name.te_$TE.blatout";

        #my $file_num         = $i + 1;
        my $te_Containing_fq = "$path/$fa_name.te_$TE.ContainingReads.fq";
        if ( -e $te_Containing_fq ) {
            $fq = $te_Containing_fq;
        }
`perl ~/bin/get_fq_of_te_trimmed_te_matching_reads.pl $path/$fa_name.te_$TE.blatout $fq $len_cutoff $mismatch_allowance > $path/$fa_name.te_$TE.flankingReads.fq `;
    }

##insert mate finder here

    my %flanking_fq;
    my @files_1 = <$path/*$mate_file_1*flankingReads.fq>;
    my @files_2        = <$path/*$mate_file_2*flankingReads.fq>;
    my @files_unpaired = <$path/*$mate_file_unpaired*flankingReads.fq>;
    if ( !@files_1 and !@files_2 and !@files_unpaired ) {
        for ( my $i = 0 ; $i < @files_1 ; $i++ ) {
            my $file_1 = $files_1[$i];
            $file_1 =~ s/$mate_file_1//;
            for ( my $j = 0 ; $j < @files_2 ; $j++ ) {
                my $file_2 = $files_2[$j];
                $file_2 =~ s/$mate_file_2//;
                if ( $file_1 eq $file_2 ) {
                    $flanking_fq{$file_1}{1} = $files_1[$i];
                    $flanking_fq{$file_1}{2} = $files_2[$j];
                    if (@files_unpaired) {
                        for ( my $k = 0 ; $k < @files_unpaired ; $k++ ) {
                            my $file_unpaired = $files_unpaired[$k];
                            if ( $file_1 eq $file_unpaired ) {
                                $flanking_fq{$file_1}{unpaired} =
                                  $files_unpaired[$k];
                                last;
                            }
                        }
                    }
                    last
                      ; #if $file_1 eq $file_2 and we are finished with unpaired go back to $i loop
                }
            }
        }
    }
    else {              ##if only unmatches files are provided
        my @files_singles = <$path/*flankingReads.fq>;
        foreach my $file ( sort @files_singles ) {
            $flanking_fq{$file}{unpaired} = $file;
        }
    }

    my @bowtie_out_files;
    my @files2merge;
    if ($mapping) {
        foreach my $key ( sort keys %flanking_fq ) {
            foreach my $type ( sort keys %{ $flanking_fq{$key} } ) {
                my $flanking_fq = $flanking_fq{$key}{$type};

                #remove and save filename part of path
                my @fq_path = split '/', $flanking_fq;
                my $fq_name = pop @fq_path;
                $fq_name =~ s/\.fq$//;
`bowtie --best -q $genome_dir/$genome_file.bowtie_build_index $flanking_fq  > $path/$target.$fq_name.bowtie.single.out `;
                push @bowtie_out_files,
                  "$path/$target.$fq_name.bowtie.single.out";
            }    #end of foreach my $type ( sort keys %{ $flanking_fq{$key} } )
            if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2} )
            {
                my $flanking_fq_1 = $flanking_fq{$key}{1};
                my $flanking_fq_2 = $flanking_fq{$key}{2};
                my @fq_path       = split '/', $flanking_fq_1;
                my $fq_name       = pop @fq_path;
                $fq_name =~ s/\.fq$//;
                if ( -s $flanking_fq_1 and -s $flanking_fq_2 ) {

                    #clean reads if both flanking.fq are non-zero file size
`~/bin/clean_pairs_memory.pl -1 $flanking_fq_1 -2 $flanking_fq_2 > $path/$fq_name.unPaired.fq`;
                }    #end of if ( -s $flanking_fq_1 and -s $flanking_fq_2 )
                if (    -s "$flanking_fq_1.matched"
                    and -s "$flanking_fq_2.matched" )
                {
`bowtie --best -q $genome_dir/$genome_file.bowtie_build_index -1 $flanking_fq_1.matched -2 $flanking_fq_2.matched > $path/$target.$fq_name.bowtie.mates.out`;
                    push @bowtie_out_files,
                      "$path/$target.$fq_name.bowtie.mates.out";
`bowtie --best -q $genome_dir/$genome_file.bowtie_build_index $path/$fq_name.unPaired.fq > $path/$target.$fq_name.bowtie.unPaired.out`;
                    push @bowtie_out_files,
                      "$path/$target.$fq_name.bowtie.unPaired.out";
                } # end of if(-s "$flanking_fq_1.matched" and -s "$flanking_fq_2.matched" )
            } #end of if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2})
        }    #end of foreach $key
        foreach my $bowtie_out (@bowtie_out_files) {

            #covert bowtie output to sam
            if ( -s $bowtie_out ) {
                `bowtie2sam.pl $bowtie_out > $bowtie_out.sam`;

                #convert sam to bam
`samtools import $genome_path.fai $bowtie_out.sam $bowtie_out.bam`;

                #sort bam
                `samtools sort $bowtie_out.bam $bowtie_out.sorted`;

                #index bam
                `samtools index $bowtie_out.sorted.bam`;
                push @files2merge, "$bowtie_out.sorted.bam";
            }    #end of if ( -s $bowtie_out )
        }    #end of foreach my $bowtie_out (@bowtie_out_files)
        my $files2merge = join " ", @files2merge;

        #merge paired and unPaired bam
        `samtools merge -f $path/$target.$TE.merged.bam $files2merge`;
`samtools sort $path/$target.$TE.merged.bam $path/$target.$TE.merged.sorted`;
        `samtools index $path/$target.$TE.merged.sorted.bam`;

        #identify mping insertion sites
`~/bin/get_TE_insertion_site.pl $path/$target.$TE.merged.sorted.bam $target $genome_path $TE $TSD{$TE} $exper`;
    }    #end of if mapping
}    #end of foreach TE fasta

print "\n####################\n#\n# output files:\n#\n####################\n\n";
my @outdirs = `ls $outdir/$top_dir`;
foreach my $out ( sort @outdirs ) {
    my @files = `ls $outdir/$top_dir/$out/*gff`;
    push @files, `ls $outdir/$top_dir/$out/*txt`;
    push @files, `ls $outdir/$top_dir/$out/*list`;
    print "##### Summary Files #####\n" if scalar @files > 0;
    print "\t## found in $outdir/$top_dir/$out ##\n";
    foreach my $file (sort @files) {
        print "$file\n";
    }
}
print "\n##### All File Descriptions #####\n";
if ($mapping) {
    print "*.te_insertion_sites.gff
 gff3 containing information about TE insertions. 
 These sites are supported by alignment of reads to both 5' and 3'
 flanking genomic sequence
        
*.te_insertion_sites.table.txt
 contains the same information about TE insertions as in the gff, 
 but in tab separtated table format.

*.te_insertion_sites.fa
 fasta containing reference genome sequence flanking the insertion site, 
 100bp-5' 100bp-3' = 100bp total for each insertion

*.te_insertion_sites.reads.list
 contains the names of reads that overlap the 5' and 3' end of TE for 
 each individual mping insertion

*.te_insertion.all.txt
 contains all possible insertion sites identified including those that
 were identified with only 5' or 3' TE flanking sequence.

*[_1|_2].te_TE.[five|three]_prime.fa
 contains the sequence matching to only TE in each read that overlaps the
 start and end of mping

*.te_TE.flankingReads[_1|_2].fq
 contains the sequence of the reads that match TE with the TE portion
 of the read removed\n";
}    #end of if mapping
else {
    print "*.te_TE.ContainingReads.fq
 contains any sequence that was found to match TE with blat\n";
}

