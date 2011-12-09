#!/usr/bin/perl -w
use strict;
use File::Spec;

my $dir = shift;

my $dir_path = File::Spec->rel2abs($dir);
$dir_path .= '/';


opendir( DIR, $dir ) || die "$!";

my %files;
foreach my $file ( readdir(DIR) ) {
    my ($volume,$directories,$filename) = File::Spec->splitpath( $file );
    
    next unless $filename =~ /fq|fastq/;
    my $fq_file = $filename;
    my $fasta_file = $filename;
    unless ($fasta_file =~ s/fq|fastq/fasta/){
   	print "$fq_file cannot be converted into a fasta file\n";
    	next;
    }

    open OUTSH, ">$fq_file.fastq2fasta.sh" or die $!;
    print OUTSH "#!/bin/bash\n\n";

    print OUTSH "fastq2fasta.pl $fq_file > $fasta_file\n";
}
