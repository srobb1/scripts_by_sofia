#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use strict;

my $genomeFasta;
my $fq_1;
my $fq_2;
my $target;
my $len_cutoff;
GetOptions(
    '1|fq_1:s'        => \$fq_1,
    '2|fq_2:s'        => \$fq_2,
    'h|help'          => \&getHelp,
);

if ( !defined $fq_1 or !defined $fq_2 ) {
    print "\n\nPlease provide 2 paired fastq files\n";
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

sub getHelp () {
    print "
usage:
./find_mping_insertionSites.pl [-g chromosome_genome_fasta][-t chromosome_to_search][-1 fastq_file_1] [-2 fastq_file_2][-h] 

options:
-1 STR          fastq file 1 (.fq or .fastq)  [no default]
-2 STR          fastq file 2 (.fq or .fastq)  [no default]
-h              this message
";

    exit 1;
}

my $genome_path = File::Spec->rel2abs($genomeFasta);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

`bowtie-build -f $genome_path bowtie_build_index`
  if !-e "bowtie_build_index.1.ebwt";
my @fq;
my @fa;
foreach my $fq ( $fq_1, $fq_2 ) {
    my $fq_path = File::Spec->rel2abs($fq);
    push @fq, $fq_path;
    if ( $fq =~ /(\/.+\/)?(\S+)\.(fq|fastq)$/ ) {
        my $fa = "$2.fa";

        #print "fa: $fa\n";
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
`perl ~/bin/get_mping_fa_for_each_read.pl $fa.blatout $fq`;
}
