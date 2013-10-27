#!/usr/bin/perl -w
use strict;
use File::Spec;


my $dir = shift;
my $dir_path = File::Spec->rel2abs($dir);



if (-d $dir_path){
  my @files = <$dir_path/*fas>;
  foreach my $file (@files){
   my @path = split /\// , $file;
   my $name = pop @path;
   my ($base) = $file =~ /(.+)\.fas$/;
   my $path = join '/' , @path;
   open SH , ">$file.sh" or die "Can't open $file.sh";
   print SH "##$file
   perl -p -e \'s/\\*\$//\' $file > $base.fa; 
   muscle -in $base.fa -out $base.aln  
   /opt/trimal/1.4/bin/trimal -automated1 -in $base.aln -out $base.trimal.aln
   perl /rhome/robb/project/Andrii/scripts/align_convert.pl $base.trimal.aln fasta $base.trimal.aln.phy phylip
";

  }
}

