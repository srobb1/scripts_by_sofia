#!/usr/bin/perl -w
use strict;
use File::Spec;
## this script uses
##use file_3 and file_4 to fill in any un_mate-ed reads in files 1 and 2

my $file_1 = shift;
my $file_2 = shift;
my $file_3 = shift;
my $file_4 = shift;
my @files_containing_te = ($file_1, $file_2);
my @files_all_records=    ($file_3, $file_4);

my %seqs;

my $count=0;
## create a hash with fq record while keeping track of input file order
foreach my $file (@files_containing_te){
	open (my $FQ_fh, "<", $file) or die "Can't open $file $!\n";
	$count++;
	print "$count: $file\n";
	while (my $fq_rec = get_FQ_record($FQ_fh)){
		my $header = get_header($fq_rec);
		my $id = $header;
		$id =~ s/\/[12]$//;
		$seqs{$id}{$count} = $fq_rec;
	}
	close $FQ_fh;
}

## if mate 1 or mate 2 is missing add it from files 3 and 4
my $this_mate = 0;
foreach my $file (@files_all_records){
	open (my $FQ_fh, "<", $file) or die "Can't open $file $!\n";
	$this_mate++;
	while (my $fq_rec = get_FQ_record($FQ_fh)){
		my $header = get_header($fq_rec);
		my $id = $header;
		$id =~ s/\/[12]$//;
		if (exists $seqs{$id}){
			my @mates = sort keys %{$seqs{$id}};
			next if scalar @mates == 2;
			next if $mates[0] == $this_mate;
			if ($mates[0] != $this_mate){
				#store it if it is the other mate
				$seqs{$id}{$this_mate}=$fq_rec;
				@mates = sort keys %{$seqs{$id}};
			}	
		}	
	}
	close $FQ_fh;
}
   
 
## print the now complete matching mates in new FQ files
## create names for the 2 new outfiles
my @outfiles;
foreach my $file (@files_containing_te){
    my $file_path = File::Spec->rel2abs( $file ) ;
    my @dirs = dir_split($file_path);
    my $filename = pop @dirs;
    my $dir_path = join '/' , @dirs;
    my @file_parts = filename_split($filename);
    push @outfiles, [ split '', $filename];
    #print "@dirs, $filename, @file_parts\n";
    #push @outfiles, "$filename.mates"; 
}
##find differnces in the two input files and delete it to make the new shuffled outfile
my $len = @{$outfiles[0]} - 1;
my @diff;
for (my $i=$len ; $i>=0  ; $i--){
        if ($outfiles[0][$i] ne $outfiles[1][$i]){
                #replace the different character and the one before it if it is a '.' or a '_'
                push @diff, $i;
                if ( $outfiles[0][$i-1] =~ /[._]/){
                        push @diff, $i-1;
                }
        }
}
foreach my $pos( sort {$b<=>$a} @diff ){
  splice @{$outfiles[0]},$pos,1;
}
my $outfile =  join ('',@{$outfiles[0]}) . ".shuffled.fq";

#print the 2 new FQ files
open (my $FQ_fh_out, ">", $outfile) or die "Can't open $outfile $!\n";
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
