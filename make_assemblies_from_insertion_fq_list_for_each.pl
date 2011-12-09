#!/usr/bin/perl -w
use strict;
use File::Spec;
## this script uses
##use file_3 and file_4 to fill in any un_mate-ed reads in files 1 and 2

my $list_file = shift;
my $dir = shift;

my @files_all_records= <$dir/*fq>;

my %seqs;
my %inserts;
## create a hash with fq ids
open (my $list_fh, '<' , $list_file) or die "Can't open $list_file $!\n";

while (my $line = <$list_fh>){
	## Chr1:174500     1:22:13443:68767:Y,1:22:13443:68767:Y
	chomp $line;
	my ($pos,$list_of_ids)= split /\t/, $line;
	#my ($chr,$coord) = split /:/,$pos; 
	$pos =~ s/:/./;
	my @ids = split ',' , $list_of_ids;
	## uniquify ids
	my %unique;
	foreach my $id (@ids){
		$id =~ s/\/[12]$//;
		$unique{$id}=1;
		$seqs{$id}{'pos'}=$pos;
	}
	my @uniq_ids = sort keys %unique;
	$inserts{$pos}=[@uniq_ids];
}
close $list_fh;
## print "looking for ", scalar (keys %seqs) , " pairs in:\n";
## if mate 1 or mate 2 is missing add it from files 3 and 4
my $this_mate = 0;
foreach my $file (sort @files_all_records){
	my $file_path = File::Spec->rel2abs( $file ) ;
	open (my $FQ_fh, "<", $file_path) or die "Can't open $file_path $!\n";
	print "\t$file_path\n";
	$this_mate++;
	my $count;
	my $rec_count;
	while (my $fq_rec = get_FQ_record($FQ_fh)){
	        $count++;
		my $header = get_header($fq_rec);
		my $id = $header;
		$id =~ s/\/[12]$//;
		$id =~ s/^@//;
		if (exists $seqs{$id}){
	         	my @mates = sort keys %{$seqs{$id}{'mate'}};
			next if scalar @mates == 2;
			#next if $mates[0] == $this_mate;
			if (!defined $mates[0] or $mates[0] != $this_mate){
		                $rec_count++;
				#store it if it is the other mate
				$seqs{$id}{'mate'}{$this_mate}=$fq_rec;
	                        print "\t\tseqs added: $rec_count\n" if $rec_count%100 == 0;	
			}
		}	
	        print "$file: seqs checked $count\n" if $count%1000000 == 0;	

	}
	close $FQ_fh;
	$this_mate = 0 if $this_mate == 2; #reset every 2 files
}
   
 
## print the now complete matching mates in new FQ files
my $file_path = File::Spec->rel2abs( $list_file ) ;
my @dirs = dir_split($file_path);
my $filename = pop @dirs;
my $dir_path = join '/' , @dirs;

open (my $sh , '>' , "velvet.sh") or die "Can't open velvet.sh for writting $!\n";
print $sh "#!/bin/bash\n\n";
foreach my $coord (keys %inserts){
  my $FQ_fh_out_2;
  my $dir = "$dir_path/$coord";
  mkdir($dir, 0755) if !-e $dir;
  my $filename = "$coord.shuffled.fq";
  open (my $FQ_fh_out, ">", "$dir/$filename") or die "Can't open $filename $!\n";
  my $unpaired = 0;
  foreach my $id ( sort @{$inserts{$coord}}){
    if (scalar keys %{$seqs{$id}{'mate'}} == 2){
     my $fq_rec_1 = join "\n", @{$seqs{$id}{'mate'}{1}};
     my $fq_rec_2 = join "\n", @{$seqs{$id}{'mate'}{2}};
    
      print $FQ_fh_out $fq_rec_1 ,"\n";
      print $FQ_fh_out $fq_rec_2 ,"\n";

    }else {
      $unpaired = 1;
      $filename = "$coord.unpaired.fq";
      open ($FQ_fh_out_2, ">", "$dir/$filename") or die "Can't open $filename $!\n";
      foreach my $mate ( keys %{$seqs{$id}{'mate'}}){
         my $fq_rec = join "\n", @{$seqs{$id}{'mate'}{$mate}};
	 print $FQ_fh_out_2 $fq_rec ,"\n";
      }
    }
 }
  close $FQ_fh_out;
  close $FQ_fh_out_2 if $unpaired;
  if ($unpaired){
    print $sh "VelvetOptimiser.pl -a -s 39 -e 61 -p $coord -f \"-fastq -shortPaired $dir/$filename -short $dir/$coord.unpaired.fq\" \n";
  }else{  
    print $sh "VelvetOptimiser.pl -a -s 39 -e 61 -p $coord -f \"-fastq -shortPaired $dir/$filename\" \n";
  }
}
close $sh;
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
