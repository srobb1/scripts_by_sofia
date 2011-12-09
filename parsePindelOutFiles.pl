#!/usr/bin/perl -w
use strict;
use File::Spec;
use Getopt::Long;

## Parses the outfile produced by Pindel 0.2.2
## *_LI *_SI *_BP *_D *_INV *_TD

## 1043    D 1     NT 0 "" ChrID chromosome12      BP 27619058     27619060        BP_range 27619058       27619063        Supports 40     + 24    - 16    S1 425  S2 1634.14      SUM_MS

my $dir = '.';
my $genomeFasta;
my $source = "pindel";

GetOptions(
    'd|dir:s' => \$dir,
    'h|help'  => \&getHelp,
);

sub getHelp () {
    print "
usage:
./pasrsePindelOutfiles.pl [-d pindel_output_directory] [-h] 

options:
-d STR		directory of Pindel Output Files  (*_D, *_SI, *_LI, *_INV, *_TD) [.]
-h 			this message
";

    exit 1;
}

my $dir_path = File::Spec->rel2abs($dir);
$dir_path .= '/';
opendir( DIR, $dir ) || die "$!";

my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

my %files;
foreach my $file ( readdir(DIR) ) {
    ##this script parses _D
    ##this script parses _SI
    ##this script parses _TD
    ##this script parses _LI

    next if $file =~ /_BP$/;
    next if $file =~ /_INV$/;

    print "$file\n";
    my ( $volume, $directories, $filename ) = File::Spec->splitpath($file);
    open INFILE, "$file" or die $!;
    while ( my $line = <INFILE> ) {
        chomp $line;
        next if $line !~ /^\d/;
        my (
            $pindel_id,
            $type_and_count_in_ref_seq_string,
            $NT_count_in_aligned_reads_and_nt_string,
            $ChrID_string,
            $BP_start_in_ref_string,
            $bp_end_in_ref_integer,
            $Supports_count_string,
            @rest
        ) = split /\t/, $line;
        if ( $file =~ /_LI$/ ) {
            my (
                $pindel_id,                  $type,
                $ChrID_string,               $bp_start_in_ref,
                $support_5prime_read_count,  $bp_end_in_ref,
                $supports_3prime_read_count, @rest
            ) = split /\t/, $line;
            my ($ref_chrID) = $ChrID_string =~ /ChrID\s+(.+)$/;
            $type = "insertion";
            my $uniqe_id =
                $ref_chrID . "_" 
              . $type . "_"
              . $bp_start_in_ref . "_"
              . $bp_end_in_ref;
            my $GFFline =
"$ref_chrID\t$source\t$type\t$bp_start_in_ref\t$bp_end_in_ref\t.\t.\t.\tID=$uniqe_id;";
	    	print "$GFFline\n"; 
        }
        else {

            ##[0]id [1]type count_in_ref_seq [2]NT count_in_aligned_seqs "nt_in_aligned_seqs" [3]ChrID name [4]BP start_in_ref [3]bp_end_in_ref [4] [5] [6] [7]
            ## 0--D 2--NT 0 ""--ChrID chromosome12--BP 227150--227153
            my ( $type, $count_in_ref ) = split /\s+/,
              $type_and_count_in_ref_seq_string;
            my ( $count_in_aligned_reads, $nt_string_in_aligned_reads ) =
              $NT_count_in_aligned_reads_and_nt_string =~ /NT\s+(.+)\s+(\S*)$/;
            my ($ref_chrID)       = $ChrID_string           =~ /ChrID\s+(.+)$/;
            my ($bp_start_in_ref) = $BP_start_in_ref_string =~ /BP\s+(\d+)$/;
            my $bp_end_in_ref     = $bp_end_in_ref_integer;
            my $supporing_reads_count =
              $Supports_count_string =~ /Supports\s+(\s+)$/;

            #print "$line\n";
            my $uniqe_id =
                $ref_chrID . "_" 
              . $type . "_"
              . $bp_start_in_ref . "_"
              . $bp_end_in_ref;
            my $column9 = "";
            if ( $type eq 'D' ) {
                $type = "deletion";
            }
            elsif ( $type eq 'SI' ) {
                $type = "insertion";
            }
            elsif ( $type eq 'TD' ) {
                my $nextLine        = <INFILE>;
                my $tandemly_dup_nt = "";
                $type = "tandem_duplication";
                if ( $nextLine =~ /([atgcn]+)/ ) {
                    $tandemly_dup_nt = $1;
                }
                $column9 = "tandemly_duplicated_nt=\"$tandemly_dup_nt\";";
            }
            if ( $type eq 'D' or $type eq 'I' ) {
                my $nextLine   = <INFILE>;
                my $deleted_nt = "";
                if ( $nextLine =~ /([atgcn]+)/ ) {
                    $deleted_nt = $1;
                }
                if ( $count_in_aligned_reads == 0 ) {    ##deletion from ref
                    $column9 = "deleted_nt=\"$deleted_nt\";";
                }
                else {    ##deletion from aligned reads => insertion into ref
                    $type    = "insertion";
                    $column9 = "inserted_nt=$nt_string_in_aligned_reads;";
                }
            }

            my $GFFline =
"$ref_chrID\t$source\t$type\t$bp_start_in_ref\t$bp_end_in_ref\t.\t.\t.\tID=$uniqe_id;";
            print "$GFFline $column9\n";
        }
    }

}

