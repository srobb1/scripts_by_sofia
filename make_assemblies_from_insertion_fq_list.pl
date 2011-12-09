#!/usr/bin/perl -w
use strict;
use File::Spec;
## this script uses
##use file_3 and file_4 to fill in any un_mate-ed reads in files 1 and 2

my $list_file = shift;
my $dir = shift;

my @files_all_records= <$dir/*fq>;

my %seqs;

## create a hash with fq ids
open (my $list_fh, '<' , $list_file) or die "Can't open $list_file $!\n";
while (my $line = <$list_fh>){
	chomp $line;
	my @ids = split ',' , $line;
	foreach my $id (@ids){
		$id =~ s/\/[12]$//;
		$seqs{$id}{0}='';
	}
}
close $list_fh;
print "looking for ", scalar (keys %seqs) , " pairs in:\n";
## if mate 1 or mate 2 is missing add it from files 3 and 4
my $this_mate = 0;
foreach my $file (sort @files_all_records){
	my $file_path = File::Spec->rel2abs( $file ) ;
	open (my $FQ_fh, "<", $file_path) or die "Can't open $file_path $!\n";
	print "\t$file_path\n";
	$this_mate++;
	while (my $fq_rec = get_FQ_record($FQ_fh)){
		my $header = get_header($fq_rec);
		my $id = $header;
		$id =~ s/\/[12]$//;
		$id =~ s/^@//;
		if (exists $seqs{$id}){
			my @mates = sort keys %{$seqs{$id}};
			next if scalar @mates > 2;
			next if $mates[0] == $this_mate;
			if ($mates[0] != $this_mate){
				#store it if it is the other mate
				$seqs{$id}{$this_mate}=$fq_rec;
			}	
		}	
	}
	close $FQ_fh;
	$this_mate = 0 if $this_mate == 2; #reset every 2 files
}
   
 
## print the now complete matching mates in new FQ files
## create names for the 2 new outfiles
    
my $file_path = File::Spec->rel2abs( $list_file ) ;
my @dirs = dir_split($file_path);
my $filename = pop @dirs;
my $dir_path = join '/' , @dirs;
$filename =~ s/\..+?$/shuffled.fq/;

#print the 2 new FQ files
open (my $FQ_fh_out, ">", $filename) or die "Can't open $filename $!\n";
foreach my $id (keys %seqs){
  my $fq_rec_1 = join "\n", @{$seqs{$id}{1}};
  my $fq_rec_2 = join "\n", @{$seqs{$id}{2}};
  print $FQ_fh_out $fq_rec_1 ,"\n";
}

#####SUBROUTINES########
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
sub get_FQ_record {
	my $file_handle = shift;
	my $ref_seq_hash = shift;
	while (my $header = <$file_handle>){
		chomp $header;
		my $seq = <$file_handle>;
		chomp $seq;
		my $qual_header = <$file_handle>;
		chomp $qual_header;
		my $qual = <$file_handle>;
		chomp $qual;
		
		die "Expected a FASTQ header line containing \'@\' but did not find it, found this instead $header\n" if substr ($header,0,1) ne '@';

		return ([$header,$seq,$qual_header,$qual]);
	}

}
sub get_header {
	my $ref2array = shift;
	return ${$ref2array}[0];
}
sub get_seq {
	my $ref2array = shift;
	return ${$ref2array}[1];
}
sub get_qual_header {
	my $ref2array = shift;
	return ${$ref2array}[2];
}
sub get_qual {
	my $ref2array = shift;
	return ${$ref2array}[3];
}
sub print_fq_record {
	my $ref2array = shift;
	return join "\n", @{$ref2array};
}
