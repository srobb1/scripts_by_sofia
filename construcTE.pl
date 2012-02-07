#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;

my $te_fq_dir   = shift or die "Please provide te containing dir";
my $all_records_dir = shift or die "Please provide all fq dir"; 
my $genome_file = shift or die "Please provide genome fasta";

my @te_fq_dir = split '/', $te_fq_dir;
my @base_dir = @te_fq_dir;
pop @base_dir;
my $base_dir        = join( '/', @base_dir );
my $bowtie_base_dir = join( '/', @base_dir );
my @fq_files        = <$te_fq_dir/*fq>;
my %seqs;
my %files;
my $construcTEr_dir = "$base_dir/construcTEr";
`mkdir -p $construcTEr_dir`;

##add all seq in TE containing fq to %seqs
foreach my $file ( sort @fq_files ) {
  my @file_name = split '/', $file;
  my $this_mate = pop @file_name;
  $this_mate =~ s/\.te_(.+)\.ContainingReads.fq$//;
  $files{$this_mate}++;
  my $te        = $1;
  my $file_path = File::Spec->rel2abs($file);
  open( my $FQ_fh, "<", $file_path ) or die "Can't open $file_path $!\n";
  while ( my $fq_rec = get_FQ_record($FQ_fh) ) {
    my $header = get_header($fq_rec);
    my $id     = $header;
    $id =~ s/\/[12]$//;
    $id =~ s/^@//;
    if ( !exists $seqs{$te}{$id}{$this_mate} ) {

      #store it if it is the other mate
      $seqs{$te}{$id}{$this_mate}{fq_rec} = $fq_rec;
      $seqs{$te}{$id}{$this_mate}{$te} = 1;
    }
  }
  close $FQ_fh;
}

##get seq of the reads and the mates that blat to the TE
foreach my $te ( keys %seqs ) {
  foreach my $file_pattern ( sort keys %files ) {
    $file_pattern =~ s/_[1|2]$//;
    my @all_records = <$all_records_dir/*$file_pattern*fq>;
    foreach my $file ( sort @all_records ) {
      my @file_name = split '/', $file;
      my $this_mate = pop @file_name;
      $this_mate =~ s/\.fq$//;
      my $file_path = File::Spec->rel2abs($file);
      open( my $FQ_fh, "<", $file_path ) or die "Can't open $file_path $!\n";
      while ( my $fq_rec = get_FQ_record($FQ_fh) ) {
        my $header = get_header($fq_rec);
        my $id     = $header;
        $id =~ s/\/[12]$//;
        $id =~ s/^@//;
        if ( exists $seqs{$te}{$id} and !exists $seqs{$te}{$id}{$this_mate} ) {
          ##store it if it is the other mate
          $seqs{$te}{$id}{$this_mate}{fq_rec} = $fq_rec;
        }
      }
      close $FQ_fh;
    }
  }

##get reads to align to genome
  my $bowtie_2_aln = "$construcTEr_dir/$te.bowtie2aln.fq";
  open BOWTIEFQ, "$bowtie_2_aln";
  foreach my $te ( keys %seqs ) {
    foreach my $id ( keys %{ $seqs{$te} } ) {
      ##next if there is only 1 mate
      next if keys %{ $seqs{$te}{$id} } == 1;
      foreach my $mate ( keys %{ $seqs{$te}{$id} } ) {
        ##2 mates
        ##for the one that originally did not align to TE
        if ( exists $seqs{$te}{$id}{$mate}{$te}
          and $seqs{$te}{$id}{$mate}{$te} != 1 )
        {
          my $seq = get_seq( $seqs{$te}{$id}{$mate}{fq_rec} );
          print BOWTIEFQ ">$mate,$id\n$seq\n";
        }
      }
    }
  }
##aln to genome
##create bowtie index
  my $bowtie_out = "$construcTEr_dir/$te.construcTEr.bowtie.out";
  if ( !-e "$genome_file.bowtie_build_index.1.ebwt" ) {
    `bowtie-build -f $genome_file $genome_file.bowtie_build_index`;
  }
`bowtie --best -q $genome_file.bowtie_build_index $bowtie_2_aln  1> $bowtie_out 2> $construcTEr_dir/bowtie.stderr`;
##parse bowtie out and record the genomic locations of alignments
  my $file_path = File::Spec->rel2abs($bowtie_out);
  open( my $BOWTIE_fh, "<", $file_path ) or die "Can't open $file_path $!\n";
  while ( my $line = <$BOWTIE_fh> ) {
    chomp $line;
    my @line = split /\t/, $line;
    my ( $name, $strand, $target, $start, $seq, $qual, $M, $mismatch ) =
      split /\t/, $line;
    my ( $this_mate, $id ) = split ',', $name;

#remove /1 or /2 from the read name
#7:12:11277:9907:Y/1     +       Chr1    22134042        TTTTTTATAAATGGATAA      DGGGGGDGGGGFGDGGGG      4       7:A>T,16:C>A
    $id =~ s/(^.+?)\/[1|2](\t.+$)/$1$2/;
    if ( exists $seqs{$te}{$id} ) {
      if ( exists $seqs{$te}{$id}{$this_mate}{$te} ) {
        delete $seqs{$te}{$id}{$this_mate}{$te};
      }
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{start}  = $start;
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{strand} = $strand;
      my @mismatches = split( ',', $mismatch );
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{mismatch} = @mismatches;
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{seq}      = $seq;
    }
  }
}
foreach my $te ( keys %seqs ) {
  my $blat_query_file = "$construcTEr_dir/$te.blat.query.fa";
  open BLATQUERY, ">$blat_query_file\n";
  foreach my $id ( keys %{ $seqs{$te} } ) {
    ##check to see if there is only 1 mate
    if ( keys %{ $seqs{$te}{$id} } == 1 ) {
      delete $seqs{$te}{$id};
      next;
    }

    #now we should have both mates
    my $te_mate     = 0;
    my $genome_mate = 0;
    foreach my $mate ( keys %{ $seqs{$te}{$id} } ) {
      foreach my $target ( keys %{ $seqs{$te}{$id}{$mate} } ) {
        if ( $target eq $te ) {
          $te_mate = $mate;
        }
        else {
          $genome_mate = 1;
        }
      }
    }
    if ( $te_mate and $genome_mate ) {
      my $seq = get_seq( $seqs{$te}{$id}{$te_mate}{$te}{fq_rec} );
      print BLATQUERY ">$te_mate,$id\n$seq\n";
    }
    else {
      delete $seqs{$te}{$id};
    }
  }
  if ( -s $blat_query_file ) {
`blat $base_dir/$te.fa $blat_query_file $construcTEr_dir/$te.construcTEr.blatout -minScore=10 -noHead`;
    open BLATOUT, "$construcTEr_dir/$te.construcTEr.blatout";
    while ( my $line = <BLATOUT> ) {
      my @line    = split /\t/, $line;
      my $strand  = $line[8];
      my $qLen    = $line[10];
      my $matches = $line[0];
      my $qName   = $line[9];
      my ( $te_mate, $id ) = split ',', $qName;
      my $tStart = $line[15] + 1;
      my $tEnd   = $line[16];
      push @{ $seqs{$te}{$id}{$te_mate}{blat_hit}{qlen} },    $qLen;
      push @{ $seqs{$te}{$id}{$te_mate}{blat_hit}{strand} },  $strand;
      push @{ $seqs{$te}{$id}{$te_mate}{blat_hit}{matches} }, $matches;
      push @{ $seqs{$te}{$id}{$te_mate}{blat_hit}{tStart} },  $tStart;
      push @{ $seqs{$te}{$id}{$te_mate}{blat_hit}{tEnd} },    $tEnd;
    }
  }
}
print Dumper \%seqs;


#####SUBROUTINES########
sub dir_split {
  my $path = shift;
  my @path = split '/', $path;
  return @path;
}

sub filename_split {
  my $file = shift;
  my @file = split /\./, $file;
  return @file;
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

sub get_FQ_record {
  my $file_handle  = shift;
  my $ref_seq_hash = shift;
  while ( my $header = <$file_handle> ) {
    chomp $header;
    my $seq = <$file_handle>;
    chomp $seq;
    my $qual_header = <$file_handle>;
    chomp $qual_header;
    my $qual = <$file_handle>;
    chomp $qual;

    die
"Expected a FASTQ header line containing \'@\' but did not find it, found this instead $header\n"
      if substr( $header, 0, 1 ) ne '@';

    return ( [ $header, $seq, $qual_header, $qual ] );
  }

}

