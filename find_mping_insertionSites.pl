#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use strict;

my $genomeFasta;
my $fq_1;
my $fq_2;
my $target;
my $len_cutoff = 10;
GetOptions(
    '1|fq_1:s'        => \$fq_1,
    '2|fq_2:s'        => \$fq_2,
    'g|genomeFasta:s' => \$genomeFasta,
    'l|len_cutoff:i'  => \$len_cutoff,
    'h|help'          => \&getHelp,
);

if ( !defined $genomeFasta ) {
    print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
    &getHelp();
}
elsif ( !defined $fq_1 or !defined $fq_2 ) {
    print "\n\nPlease provide 2 paired fastq files\n";
    &getHelp();
}
unless ( -e $genomeFasta ) {
    print "$genomeFasta does not exist. Check file name.\n";
    &getHelp();
}
unless ( -e $fq_1 ) {
    print "$fq_1 does not exist. Check file name.\n";
    &getHelp();
}
unless ( -e $fq_2 ) {
    print "$fq_2 does not exist. Check file name.\n";
    &getHelp();
}
if ( !defined $target ) {
	open FASTA, $genomeFasta;
	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ />(\S+)/){
			$target = $1;
			last;
		}
	}
}
sub getHelp () {
    print "
usage:
./find_mping_insertionSites.pl [-g chromosome_genome_fasta][-1 fastq_file_1] [-2 fastq_file_2][-h] 

options:
-g STR          single chromosome genome fasta file path [no default]
-1 STR          fastq file 1 (.fq or .fastq)  [no default]
-2 STR          fastq file 2 (.fq or .fastq)  [no default]
-l INT		len cutoff for the mping trimmed reads to be aligned [10] 
-h              this message
";

    exit 1;
}

my $genome_path = File::Spec->rel2abs($genomeFasta);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

`bowtie-build -f $genome_path $target.bowtie_build_index`
  if !-e "$target.bowtie_build_index.1.ebwt";
my @fq;
my @fa;
foreach my $fq ( $fq_1, $fq_2 ) {
    my $fq_path = File::Spec->rel2abs($fq);
    push @fq, $fq_path;
    if ( $fq =~ /(\/.+\/)?(\S+)\.(fq|fastq)$/ ) {
        my $fa = "$2.fa";

        push @fa, $fa;
        if ( !-e $fa ) {
            open INFQ,  $fq_path or die $1;
            open OUTFA, ">$fa"   or die $1;

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

#create mping.fa
if ( !-e "mping.fa" ) {
    open MPING, ">mping.fa" or die $!;
    print MPING
">mping gi|22830894|dbj|AB087615.1| Oryza sativa Japonica Group transposon mPing/miniSNOOPY, complete sequence
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATG
ATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTT
TCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGT
CCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAA
CTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGT
TTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATT
GTGACTGGCC";
    close MPING;
}

#blat fa files against mping.fa
my @flanking_fq;
for ( my $i = 0 ; $i < 2 ; $i++ ) {
    my $fa = $fa[$i];
    my $fq = $fq[$i];
    `blat -minScore=10 -tileSize=7 mping.fa $fa $fa.blatout`
      if !-e "$fa.blatout";
    my $file_num = $i + 1;
    $fa =~ /(\S+)\.fa/;	
    my $mpingContaining_fq = "$1.mpingContainingReads.fq";
    if (-e $mpingContaining_fq){
	$fq = $mpingContaining_fq;
    }
`perl ~/bin/get_fq_of_mpingTrimmed_mpingMatching_reads.pl $fa.blatout $fq $len_cutoff > mpingFlankingReads_$file_num.fq`;
    push @flanking_fq, "mpingFlankingReads_$file_num.fq";
}

#clean reads
`~/bin/clean_pairs_memory.pl -1 $flanking_fq[0] -2 $flanking_fq[1] > mpingFlanking_unPaired.fq`;

#align mpingFlanking reads to genome fasta
`bowtie --best -q $target.bowtie_build_index -1 $flanking_fq[0].matched -2 $flanking_fq[1].matched > bowtie.out`;
`bowtie --best -q $target.bowtie_build_index mpingFlanking_unPaired.fq > bowtie_unPaired.out`;

#create an index of genome fasta
`samtools faidx $genome_path`;

#covert bowtie output to sam
`bowtie2sam.pl bowtie.out > mpingFlankingReads.bowtie.aln.sam`;
`bowtie2sam.pl bowtie_unPaired.out > mpingFlankingReads.bowtie.aln.unPaired.sam`;

#convert sam to bam
`samtools import $genome_path.fai mpingFlankingReads.bowtie.aln.sam mpingFlankingReads.bowtie.aln.bam`;
`samtools import $genome_path.fai mpingFlankingReads.bowtie.aln.unPaired.sam mpingFlankingReads.bowtie.aln.unPaired.bam`;

#sort bam
`samtools sort mpingFlankingReads.bowtie.aln.bam mpingFlankingReads.bowtie.aln.sorted`;
`samtools sort mpingFlankingReads.bowtie.aln.unPaired.bam mpingFlankingReads.bowtie.aln.unPaired.sorted`;

#index bam
`samtools index mpingFlankingReads.bowtie.aln.sorted.bam`;
`samtools index mpingFlankingReads.bowtie.aln.unPaired.sorted.bam`;

#merge paired and unPaired bam
`samtools merge -f mpingFlankingReads.bowtie.aln.merged.bam mpingFlankingReads.bowtie.aln.sorted.bam mpingFlankingReads.bowtie.aln.unPaired.sorted.bam`;
`samtools sort mpingFlankingReads.bowtie.aln.merged.bam mpingFlankingReads.bowtie.aln.merged.sorted`;
`samtools index mpingFlankingReads.bowtie.aln.merged.sorted.bam`;

#identify mping insertion sites
`~/bin/get_mping_insertion_site.pl mpingFlankingReads.bowtie.aln.merged.sorted.bam $target $genome_path`;

print "\n\noutput files:\n\n";
print
"$target.mping_insertion_sites.gff\tgff3 containing information about mping insertions. These sites are supposrted by alignment of reads to both 5' and 3' flanking genomic sequence\n";
print
"$target.mping_insertion_sites.table.txt\tcontains the same information about mping insertions as in the gff, but in tab separtated table format.\n";
print
"$target.mping_insertion_sites.fa\tfasta containing reference genome sequence flanking the insertion site, 50bp-5' 50bp-3' = 100bp total for each insertion\n";
print
"$target.mping_insertion_sites.reads.list\tcontains the names of reads that overlap the 5' and 3' end of mping for each individual mping insertion\n";
print $target
  . "[_1|_2].mping_[five|three]_prime.fa\tcontains the sequence matching to only mping in each read that overlaps the start and end of mping\n";
print $target
  . "[_1|_2].mpingContainingReads.fq\tcontains any sequence that was found to match mping with blat\n";
print
"$target.mping_insertion_sites.all.txt\tcontains all possible insertion sites identified including those that were identified with only 5' or 3' mping flanking sequence.\n";
print
"mpingFlankingReads[_1|_2].fq\tcontains the sequence of the reads that match mping with the mping portion of the read removed\n";
