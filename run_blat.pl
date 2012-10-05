#!/usr/bin/perl -w
use strict;
use File::Spec;
my $fa_dir = shift;
my $db=shift;
my $blatout_dir = shift;

my $cwd = `pwd`;
$cwd =~ s/\s//g;



if (!defined $blatout_dir){
  $blatout_dir = "$cwd/blatout";
}else{
  my @blatout_dir = split '/' , $blatout_dir;
  my $top = pop @blatout_dir;
  $blatout_dir = "$cwd/$top";
}
if (!-e $blatout_dir){
  `mkdir -p $blatout_dir`;
}

my $full_fa_dir = File::Spec->rel2abs($fa_dir);
my $db_path =  File::Spec->rel2abs($db);
die "Please provide a fasta to be used as a database\n" if !defined $db_path or !-e $db_path;
 
my @fa_files = <$full_fa_dir/*fa>;

foreach my $fa_file (@fa_files){
  chomp $fa_file;
  my @path = split /\// , $fa_file;
  my $file_name = pop @path;
  my $path = join '/' , @path;

  open OUTSH, ">$cwd/$file_name.blat.sh" or die "Can't open $cwd/$file_name.blat.sh for writting\n";
  print OUTSH "tmp=\`mktemp -d -p /scratch\`\n";
  print OUTSH "cd \$tmp\n"; 
  print OUTSH "blat $db_path $fa_file $file_name.vs.uniq.blatout\n"; 
  print OUTSH "mv $file_name.vs.uniq.blatout $blatout_dir/$file_name.vs.uniq.blatout\n"; 
  print OUTSH "rm -rf \$tmp\n"; 
}
