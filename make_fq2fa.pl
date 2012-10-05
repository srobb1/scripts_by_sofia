#!/usr/bin/perl -w
use strict;
use File::Spec;
my $fq_dir = shift;
my $full_fq_dir = File::Spec->rel2abs($fq_dir);


opendir( DIR, $full_fq_dir ) || die "$!";
foreach my $fq_file ( readdir(DIR) ) {
  next unless $fq_file =~ /\.(fastq|fq)$/;
  my $ext = $1;
  my $fa_file = $fq_file;
  $fa_file =~ s/$ext/fa/;
  my @path = split '/' , $fa_file;  
  my $pre = pop @path;
  $pre =~ s/\.fa//; 
  open OUTSH, ">$pre.fq2fa.sh";
  print OUTSH "tmp=\`mktemp -d -p /scratch\`
cd \$tmp
ln -s $full_fq_dir/$fq_file tmp.fastq 
/home_stajichlab/robb/src/RelocaTE_pages/scripts/relocaTE_fq2fa.pl tmp.fastq fasta
mv fasta $full_fq_dir/$fa_file
rm -rf \$tmp " if !-e $fa_file;
}
