#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use Bio::SeqIO;

use strict;
if ( !defined @ARGV ) {
    &getHelp();
}

my $te_fasta;
my $len_cutoff         = 10;
my $mismatch_allowance = 0;
my $fa_dir;

my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

GetOptions(
    'd|fa_dir:s'      => \$fa_dir,
    't|te_fasta:s'    => \$te_fasta,
    'l|len_cutoff:i'  => \$len_cutoff,
    'm|mismatch:f'    => \$mismatch_allowance,
    'h|help'          => \&getHelp,
);

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
my @fa_files;
my %fa_files;
if ( !defined $fa_dir ) {
    print "\n\nPlease provide a directory of fasta files\n";
    &getHelp();
}
elsif ( !-d $fa_dir ) {
    print
"\n\nCheck the spelling or location of $fa_dir, Please provide a directory of fasta files\n";
    &getHelp();
}
else {

    my $fa_path = File::Spec->rel2abs($fa_dir);
    #print "fa_path: $fa_path\n";
    @fa_files = <$fa_path/*fa>;
    my @fasta_files = <$fa_path/*fasta>;

    push @fa_files, @fasta_files;

    foreach my $file (@fa_files){
	#print "fa: $file\n";
    }
    #checking for more than 0 files
    if ( scalar @fa_files == 0 ) {
        print
"Must provide at least 1 fasta file\n";
        &getHelp();
    }
}

sub getHelp {
    print "
usage:
./find_TE_insertionSites.pl [-t TE_fasta_file][-d dir_of_fq][-m mismatch_allowance][-h] 

options:
-t STR          fasta containing 1 or more nucleotide sequences of transposable elements with TSD in the desc [no default]
-d STR          directory of paired fastq files (paired _1.fq & _2.fq) (.fq or .fastq is acceptable)  [no default]
-l INT		len cutoff for the te trimmed reads to be aligned [10] 
-m FRACTION	mismatch allowance for alignment to TE (int, ex 0.1) [0] 
-h              this message

SAMPLE TE FASTA
>mping	TTA
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

my $te_path = File::Spec->rel2abs($te_fasta);

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
	#my $te_dir = "$current_dir/$id";
        push @te_fastas, "$current_dir/$te_file";
	#`mkdir -p $te_dir`;
        open( OUTFASTA, ">$current_dir/$te_file" ) or die "$!\n";
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
    my @path = split '/' , $te_path;
    my $te_fasta = pop @path;
    my $path = join '/',@path;
    my $TE = $te_fasta;
    $TE =~ s/\.fa//;

    #blat fa files against te.fa
    my @flanking_fa;
    my $fa_file_count = scalar @fa_files;
    for ( my $i = 0 ; $i < $fa_file_count ; $i++ ) {
        my $fa = $fa_files[$i];
	#print "blat-ing $fa against  $path/$te_fasta\n";
	#remove and save filename part of path
        my @fa_path = split '/' , $fa;
	my $fa_name = pop @fa_path;
	$fa_name =~ s/\.fa$//;
	
        `blat -minScore=10 -tileSize=7 $path/$te_fasta $fa $path/$fa_name.te_$TE.blatout`
          if !-e "$path/$fa_name.te_$TE.blatout";
        my $file_num         = $i + 1;
        my $te_Containing_fq = "$path/$fa_name.te_$TE.ContainingReads.fq";
        #if ( -e $te_Containing_fq ) {
        #    $fq = $te_Containing_fq;
        #}
`perl ~/bin/get_fa_of_te_trimmed_te_matching_reads.pl $path/$fa_name.te_$TE.blatout $fa $len_cutoff $mismatch_allowance > $path/$fa_name.te_$TE.flankingReads.fa `;
        push @flanking_fa, "$path/$fa_name.te_$TE.flankingReads.fa";
    }

foreach my $file (@flanking_fa){
	my @dir = dir_split ($file);
	my $filename = pop @dir;
        my $branch = pop @dir;
	$filename =~ s/\.(fa|fasta)//;
	my $seqIO_obj = Bio::SeqIO->new(-file => $file , -format=>'fasta');
        while (my $seq_obj = $seqIO_obj->next_seq()){
	        my $seq = $seq_obj->seq;
                my $id = $seq_obj->id;
                my $desc = $seq_obj->desc;
                #print "found $id in $file\n";
=cut ##think that i do not want to revcomp here, i want to compare the seq from the TTA/TAA to the left and to the right
	 	if ($seq =~ /(TTA|TAA)$/){
			#revcomp
			$seq = $seq_obj->revcom->seq;
                        
		}
=cut
		print ">$branch.$filename.$TE.$id $desc\n$seq\n";
	}
}
}
sub dir_split{
        my $path = shift;
        my @path = split '/', $path;
        return @path;
}
sub filename_split {
        my $file = shift;
        my @file = split /\./ , $file;
        return @file;
}
