#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;

my $all_records_dir = shift;
my $genome_file     = shift;
my $te_fasta        = shift;
my $working_dir     = shift;
my $regex_file      = shift;
my %seqs;

##get the regelar expression patterns for mates and for the TE
##when passed on the command line as an argument, even in single
##quotes I lose special regex characters
open INREGEX, "$regex_file" or die "$!";
my $mate_1_pattern;
my $mate_2_pattern;

while ( my $line = <INREGEX> ) {
  chomp $line;
  ( $mate_1_pattern, $mate_2_pattern ) = split /\t/, $line;
  print "$mate_1_pattern--$mate_2_pattern\n";
}

## store ranges in TE fasta in %ranges
print "Storing Ranges\n";
my %ranges;
open TE_FA, "$te_fasta" or die "Can't open $te_fasta\n";

while ( my $line = <TE_FA> ) {
  chomp $line;
  if ( $line =~ /^>(\S+)/ ) {
    my $te = $1;
    while ( $line =~ /range=(\S+):(\d+)..(\d+)/ig ) {
      my $range_name = $1;
      my $s          = $2;
      my $e          = $3;

      $ranges{$te}{$range_name} = range_create( $s, $e );
    }
  }
}
print "Finding Blat files in $working_dir/blat_output\n";
my @blat_files = <$working_dir/blat_output/*blatout>;
foreach my $blat_file (@blat_files) {
  ##blat parser
  my @file_path = split '/', $blat_file;
  my $file_name = pop @file_path;
  my $FA        = $file_name;
  $FA =~ s/te_(.+).blatout/fa/;
  my $te      = $1;
  my $te_mate = $FA;
  $te_mate =~ s/\.fa//;

  open INBLAT, $blat_file, or die "Please provide a blat output file\n";

  <INBLAT>;    #get rid of blat header lines
  <INBLAT>;
  <INBLAT>;
  <INBLAT>;
  <INBLAT>;
  my %seq_storage;
  while ( my $line = <INBLAT> ) {
    my @line        = split /\t/, $line;
    my $matches     = $line[0];
    my $mismatches  = $line[1];
    my $qBaseInsert = $line[5];
    my $tBaseInsert = $line[7];
    my $strand      = $line[8];
    my $qName       = $line[9];
    my $qLen        = $line[10];
    my $tLen        = $line[14];
    my $tStart      = $line[15] + 1;
    my $tEnd        = $line[16];
    my $id          = $qName;
    my $aln_bp      = $matches + $qBaseInsert + $mismatches;
    ## if the hit overlaps the edge of the TE and not most of the read is aligning
    ## throw it out
    next
      if (( $tStart == 1 or $tEnd == $tLen )
      and ($aln_bp) <= ( $qLen * .98 ) );
    my $add = 0;
    if ( exists $seqs{$te}{$id}{$te_mate}{blat_hit}{matches} ) {
      my $stored_matches = $seqs{$te}{$id}{$te_mate}{blat_hit}{matches};
      my $stored_mm      = $seqs{$te}{$id}{$te_mate}{blat_hit}{mismatches};
      my $stored_qBI     = $seqs{$te}{$id}{$te_mate}{blat_hit}{qBaseInsert};
      my $stored_aln_bp  = $stored_matches + $stored_qBI + $stored_mm;
      ## if blat hit for this read already exists, is it better?
      if ( ($aln_bp) > ($stored_aln_bp) ) {
        $add = 1;
        delete $seqs{$te}{$id}{$te_mate}{blat_hit};
      }
      else {
        ##don't do anything, don't add new one, don't delete old one
      }
    }
    else {
      ##doesnt exist so add it
      $add = 1;
    }
    if ($add) {
      $seqs{$te}{$id}{$te_mate}{blat_hit}{qlen}        = $qLen;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{qBaseInsert} = $qBaseInsert;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{tBaseInsert} = $tBaseInsert;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{strand}      = $strand;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{matches}     = $matches;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{mismatches}  = $mismatches;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{tStart}      = $tStart;
      $seqs{$te}{$id}{$te_mate}{blat_hit}{tEnd}        = $tEnd;
      $seqs{$te}{$id}{$te_mate}{$te}                   = 1;
      my $other_mate = $te_mate;

      if ( $other_mate =~ /$mate_1_pattern/ ) {
        $other_mate =~ s/$mate_1_pattern/$mate_2_pattern/;
      }
      else {
        $other_mate =~ s/$mate_2_pattern/$mate_1_pattern/;
      }
      $seq_storage{$te_mate} .= "$id,";
      $seq_storage{$other_mate} .= "$id,";
    }
  }
  foreach my $mate ( keys %seq_storage ) {
    my $fastacmd_str = $seq_storage{$mate};
    ##make this a subroutine
    #$fastacmd_str =~ s/,$//;
    my $seq_recs =
      `fastacmd -d $all_records_dir/$mate.fa -s $fastacmd_str`;
    if ( defined $seq_recs ) {
      my @seq_recs = split />/, $seq_recs;
      ## get rid of first empty record
      shift @seq_recs; 
      ## each record = id\nseq\nseq\n
      foreach my $rec (@seq_recs) {
        my ( $header, @seq ) = split /\n/, $rec;
        my ($id) = split /\s+/, $header;
        $id =~ s/lcl\|//;
        my $seq = join '', @seq;
        $seqs{$te}{$id}{$te_mate}{seq} = $seq;
      }
    }
  }
}
print Dumper \%seqs;
##checking to see if these files already exist
#my $out_fa_path = "$working_dir/te_containing_fa";
#if (! -e "$out_fa_path/$FA.te_$TE.ContainingReads.fa" ) {
#  open OUTFQ,  ">$out_fq_path/$FA.te_$TE.ContainingReads.fa";
#  my $fastacmd_str = '';
#  foreach my $id (keys %{$seqs{$TE}}){
#    $fastacmd_str .= "$id,";
#  }
#  print OUTFQ `fastacmd -d $all_records_dir/$FA.fa -p F -s $fastacmd_str`;
#  `formatdb -i $out_fq_path/$FA.te_$TE.ContainingReads.fa -p F -o T`;
#}

#  $relocaTE =~ s/\/$//;
#  $te_fq_dir = "$relocaTE/te_containing_fq";

#my @te_fq_dir = split '/', $te_fq_dir;
#my @base_dir = @te_fq_dir;
#pop @base_dir;
#my $base_dir        = join( '/', @base_dir );
#my $bowtie_base_dir = join( '/', @base_dir );
#my @fq_files        = <$te_fq_dir/*fq>;
#my %seqs;
#my %files;
#my $construcTEr_dir = "$base_dir/construcTEr";
#`mkdir -p $construcTEr_dir`;

##add all seq in TE containing fq to %seqs
#foreach my $file ( sort @fq_files ) {
#  my @file_name = split '/', $file;
#  my $this_mate = pop @file_name;
#  $this_mate =~ s/\.te_(.+)\.ContainingReads.fq$//;
#  my $te           = $1;
#  my $file_pattern = $this_mate;
#  $file_pattern =~ s/_[1|2]$//;
#  $files{$file_pattern}++;
#  my $file_path = File::Spec->rel2abs($file);
#  open( my $FQ_fh, "<", $file_path ) or die "Can't open $file_path $!\n";

#  while ( my $fq_rec = get_FQ_record($FQ_fh) ) {
#    my $header = get_header($fq_rec);
#    my $id     = $header;
#    $id =~ s/\/[12]$//;
#    $id =~ s/^@//;
#    if ( !exists $seqs{$te}{$id}{$this_mate} ) {

#store it if it is the other mate
#     $seqs{$te}{$id}{$this_mate}{fq_rec} = $fq_rec;
#     $seqs{$te}{$id}{$this_mate}{$te} = 1;
#   }
# }
# close $FQ_fh;
#}

##get seq of the reads and the mates that blat to the TE
#foreach my $file_pattern ( sort keys %files ) {
#  my @all_records = <$all_records_dir/*$file_pattern*fq>;
#  foreach my $file ( sort @all_records ) {
#    my @file_name    = split '/', $file;
#    my $this_mate    = pop @file_name;
#    my $te_and_mates = "$construcTEr_dir/$this_mate";
#    $this_mate =~ s/\.fq$//;
#    my $file_path;
#    my $cherry_picked = 0;
#    if ( -e $te_and_mates and -s $te_and_mates ) {
#      $file_path     = File::Spec->rel2abs($te_and_mates);
#      $cherry_picked = 1;
#    }
#    else {
#      $file_path = File::Spec->rel2abs($file);
#      open OUTFQ, ">$te_and_mates"
#        or die "Can't open $te_and_mates for writing , $!\n";
#    }
#    open( my $FQ_fh, "<", $file_path ) or die "Can't open $file_path $!\n";
#    while ( my $fq_rec = get_FQ_record($FQ_fh) ) {
#      my $header = get_header($fq_rec);
#      my $id     = $header;
#      $id =~ s/\/[12]$//;
#      $id =~ s/^@//;
#      foreach my $te ( keys %seqs ) {
#        if ( exists $seqs{$te}{$id} and !exists $seqs{$te}{$id}{$this_mate} ) {
#          ##store it if it is the other mate
#          $seqs{$te}{$id}{$this_mate}{fq_rec} = $fq_rec;
#          print OUTFQ print_fq_record($fq_rec), "\n" if !$cherry_picked;
#        }
#      }
#    }
#    close OUTFQ if !$cherry_picked;
#    close $FQ_fh;
#  }
#}
#
##get reads to align to genome
foreach my $te ( keys %seqs ) {

  #my $bowtie_2_aln = "$construcTEr_dir/$te.construcTEr.bowtie2aln.fq";
  #open BOWTIEFQ, ">$bowtie_2_aln" or die "Can't open $bowtie_2_aln, $!\n";
  my $bowtie_2_aln = "$working_dir/$te.construcTEr.bowtie2aln.fa";
  open BOWTIEFA, ">$bowtie_2_aln" or die "Can't open $bowtie_2_aln, $!\n";
  foreach my $id ( keys %{ $seqs{$te} } ) {
    ##next if there is only 1 mate
    next if keys %{ $seqs{$te}{$id} } == 1;
    foreach my $mate ( keys %{ $seqs{$te}{$id} } ) {
      ##2 mates
      ##for the one that originally did not align to TE
      if ( exists $seqs{$te}{$id}{$mate}
        and !exists $seqs{$te}{$id}{$mate}{$te} )
      {
        my $seq = $seqs{$te}{$id}{$mate}{seq};
        print BOWTIEFA "$mate,$id\n$seq\n";

        #my $seq         = get_seq( $seqs{$te}{$id}{$mate}{fq_rec} );
        #my $qual        = get_qual( $seqs{$te}{$id}{$mate}{fq_rec} );
        #my $qual_header = get_qual_header( $seqs{$te}{$id}{$mate}{fq_rec} );
        #print BOWTIEFQ '@', "$mate,$id\n$seq\t$qual_header\n$qual\n";
      }
    }
  }
##aln to genome
##create bowtie index
  my $bowtie_out = "$working_dir/$te.construcTEr.bowtie.out";
  if ( !-e "$genome_file.bowtie_build_index.1.ebwt" ) {
    `bowtie-build -f $genome_file $genome_file.bowtie_build_index`;
  }
`bowtie --best -a -v 2 -q $genome_file.bowtie_build_index $bowtie_2_aln  1> $bowtie_out 2> $working_dir/bowtie.stderr`;
##parse bowtie out and record the genomic locations of alignments
  my $file_path = File::Spec->rel2abs($bowtie_out);
  open( my $BOWTIE_fh, "<", $file_path ) or die "Can't open $file_path $!\n";
  while ( my $line = <$BOWTIE_fh> ) {
    chomp $line;
    my @line = split /\t/, $line;
    my ( $name, $strand, $target, $start, $seq, $qual, $M, $mismatch ) =
      split /\t/, $line;
    my ( $this_mate, $id ) = split ',', $name;
    ##remove /1 or /2 from the read name
    ##7:12:11277:9907:Y/1     +       Chr1    22134042        TTTTTTATAAATGGATAA      DGGGGGDGGGGFGDGGGG      4       7:A>T,16:C>A
    $id =~ s/(^.+?)\/[1|2](\t.+$)/$1$2/;
    if ( exists $seqs{$te}{$id} ) {
      ##if an aln to genome exsts, delete record that says it aligns to TE
      if ( exists $seqs{$te}{$id}{$this_mate}{$te} ) {
        delete $seqs{$te}{$id}{$this_mate}{$te};
      }
      my $end = $start + ( length $seq ) - 1;
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{"$start..$end($mismatch)"}
        {start} = $start;
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{"$start..$end($mismatch)"}
        {strand} = $strand;
      my @mismatches = split( ',', $mismatch );
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{"$start..$end($mismatch)"}
        {mismatch} = @mismatches;
      $seqs{$te}{$id}{$this_mate}{aln}{$target}{"$start..$end($mismatch)"}
        {seq} = $seq;
    }
  }
}

#my %ranges;
#foreach my $te ( keys %seqs ) {
#  open TE_FA, "$base_dir/$te.fa" or die "Can't open $base_dir/$te.fa\n";
#  while (my $line = <TE_FA>){
#    chomp $line;
#    if ($line =~ /^>(\S+)/){
#      my $te = $1;
#      while ($line =~ /range=(\S+):(\d+)..(\d+)/ig){
#        my $range_name = $1;
#        my $s = $2;
#        my $e = $3;
#
#       $ranges{$te}{$range_name} = range_create($s, $e);
#      }
#    }
#  }
#  my $blat_query_file = "$construcTEr_dir/$te.blat.query.fa";
#  open BLATQUERY, ">$blat_query_file\n";
#  foreach my $id ( keys %{ $seqs{$te} } ) {
#    ##check to see if there is only 1 mate
#    if ( keys %{ $seqs{$te}{$id} } == 1 ) {
#      delete $seqs{$te}{$id};
#      next;
#    }
#    ##now we should have both mates
#    my $te_mate     = 0;
#    my $genome_mate = 0;
#    foreach my $mate ( keys %{ $seqs{$te}{$id} } ) {
#      foreach my $target ( keys %{ $seqs{$te}{$id}{$mate} } ) {
#        if ( $target eq $te ) {
#          $te_mate = $mate;
#        }
#        else {
#          $genome_mate = 1;
#        }
#      }
#    }
#    if ( $te_mate and $genome_mate ) {
#      my $seq = get_seq( $seqs{$te}{$id}{$te_mate}{fq_rec} );
#      print BLATQUERY ">$te_mate,$id\n$seq\n";
#    }
#    else {
#      delete $seqs{$te}{$id};
#    }
#  }
#  if ( -s $blat_query_file ) {
#`blat $base_dir/$te.fa $blat_query_file $construcTEr_dir/$te.construcTEr.blatout -minScore=10 -noHead`;
#    open BLATOUT, "$construcTEr_dir/$te.construcTEr.blatout";
#    while ( my $line = <BLATOUT> ) {
#      my @line        = split /\t/, $line;
#      my $matches     = $line[0];
#      my $mismatches  = $line[1];
#      my $qBaseInsert = $line[5];
#      my $tBaseInsert = $line[7];
#      my $strand      = $line[8];
#      my $qName       = $line[9];
#      my $qLen        = $line[10];
#      my $tLen        = $line[14];
#      my $tStart      = $line[15] + 1;
#      my $tEnd        = $line[16];
#      my ( $te_mate, $id ) = split ',', $qName;
#      my $add = 0;
#
#      my $aln_bp         = $matches + $qBaseInsert + $mismatches ;
#      ## if the hit overlaps the edge of the TE and not most of the read is aligning
#      ## throw it out
#      next if (( $tStart == 1 or $tEnd == $tLen ) and ( $aln_bp ) <= ( $qLen * .98 ) );
#      if ( exists $seqs{$te}{$id}{$te_mate}{blat_hit}{matches} ) {
#        my $stored_matches = $seqs{$te}{$id}{$te_mate}{blat_hit}{matches};
#        my $stored_mm      = $seqs{$te}{$id}{$te_mate}{blat_hit}{mismatches};
#        my $stored_qBI     = $seqs{$te}{$id}{$te_mate}{blat_hit}{qBaseInsert};
#        my $stored_aln_bp  = $stored_matches + $stored_qBI + $stored_mm ;
#        ## if blat hit for this read already exists, is it better?
#        if ( ( $aln_bp ) > ( $stored_aln_bp ) ) {
#          $add = 1;
#          delete $seqs{$te}{$id}{$te_mate}{blat_hit};
#        }
#        else {
#          ##don't do anything, don't add new one, don't delete old one
#        }
#      }
#      else {
#        ##doesnt exist so add it
#        $add = 1;
#      }
#      if ($add) {
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{qlen}        = $qLen;
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{qBaseInsert} = $qBaseInsert;
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{tBaseInsert} = $tBaseInsert;
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{strand}      = $strand;
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{matches}     = $matches;
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{mismatches}  = $mismatches;
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{tStart}      = $tStart;
#        $seqs{$te}{$id}{$te_mate}{blat_hit}{tEnd}        = $tEnd;
#      }
#    }
#  }
#}
#
my %construcTEr;
foreach my $te ( keys %seqs ) {
  foreach my $id ( keys %{ $seqs{$te} } ) {
    my $te_hit     = 0;
    my $genome_aln = 0;
    foreach my $mate ( keys %{ $seqs{$te}{$id} } ) {
      if ( exists $seqs{$te}{$id}{$mate}{blat_hit} ) {
        ##my $matches         = $seqs{$te}{$id}{$mate}{blat_hit}{matches};
        ##my $q_len           = $seqs{$te}{$id}{$mate}{blat_hit}{qlen};
        ##my $percent_aligned = $matches / $q_len;
        ##if ( $percent_aligned >= .95 ) {
        $te_hit = 1;
        ##}
      }
      if ( exists $seqs{$te}{$id}{$mate}{aln} ) {
        foreach my $target ( sort keys %{ $seqs{$te}{$id}{$mate}{aln} } ) {
          foreach
            my $id_str ( sort keys %{ $seqs{$te}{$id}{$mate}{aln}{$target} } )
          {
            my $mismatches =
              $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{mismatch};
            my $seq_len =
              length $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{seq};
            my $percent_mm = $mismatches / $seq_len;
            if ( $percent_mm < 0.1 ) {
              $genome_aln = 1;
            }
            else {
              delete $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str};
            }
          }
        }
      }
    }
    if ( $te_hit and $genome_aln ) {
      my %storage;
      foreach my $mate ( sort keys %{ $seqs{$te}{$id} } ) {
        my $fq_obj = $seqs{$te}{$id}{$mate}{fq_rec};
        my $seq    = get_seq($fq_obj);
        if ( exists $seqs{$te}{$id}{$mate}{blat_hit} ) {
          my $tStart     = $seqs{$te}{$id}{$mate}{blat_hit}{tStart};
          my $qlen       = $seqs{$te}{$id}{$mate}{blat_hit}{qlen};
          my $matches    = $seqs{$te}{$id}{$mate}{blat_hit}{matches};
          my $mismatches = $seqs{$te}{$id}{$mate}{blat_hit}{mismatches};
          my $tEnd       = $seqs{$te}{$id}{$mate}{blat_hit}{tEnd};
          my $strand     = $seqs{$te}{$id}{$mate}{blat_hit}{strand};
          $storage{TE}{$te}{start}      = $tStart;
          $storage{TE}{$te}{end}        = $tEnd;
          $storage{TE}{$te}{mismatches} = $mismatches;
          $storage{TE}{$te}{seq}        = $seq;
          $storage{TE}{$te}{strand}     = $strand;
        }
        elsif ( exists $seqs{$te}{$id}{$mate}{aln} ) {
          foreach my $target ( sort keys %{ $seqs{$te}{$id}{$mate}{aln} } ) {
            foreach
              my $id_str ( sort keys %{ $seqs{$te}{$id}{$mate}{aln}{$target} } )
            {
              my $start = $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{start};
              my $strand =
                $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{strand};
              my $mismatches =
                $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{mismatch};
              my $seq_len =
                length $seqs{$te}{$id}{$mate}{aln}{$target}{$id_str}{seq};
              my $end = $start + $seq_len - 1;
              $storage{genome}{$target}{start}      = $start;
              $storage{genome}{$target}{end}        = $end;
              $storage{genome}{$target}{mismatches} = $mismatches;
              $storage{genome}{$target}{strand}     = $strand;
              $storage{genome}{$target}{seq}        = $seq;
            }
          }
        }
      }
      foreach my $te ( sort keys %{ $storage{TE} } ) {
        my $te_start      = $storage{TE}{$te}{start};
        my $te_read_seq   = $storage{TE}{$te}{seq};
        my $te_end        = $storage{TE}{$te}{end};
        my $te_strand     = $storage{TE}{$te}{strand};
        my $te_mismatches = $storage{TE}{$te}{mismatches};
        foreach my $target ( sort keys %{ $storage{genome} } ) {
          my $start           = $storage{genome}{$target}{start};
          my $genome_read_seq = $storage{genome}{$target}{seq};
          my $end             = $storage{genome}{$target}{end};
          my $strand          = $storage{genome}{$target}{strand};
          my $mismatches      = $storage{genome}{$target}{mismatches};
          $construcTEr{$target}{"$start..$end"}{$target}{start} = $start;
          $construcTEr{$target}{"$start..$end"}{$target}{end}   = $end;
          $construcTEr{$target}{"$start..$end"}{$target}{mismatch} =
            $mismatches;
          $construcTEr{$target}{"$start..$end"}{$target}{strand} = $strand;
          $construcTEr{$target}{"$start..$end"}{$target}{id}     = $id;
          $construcTEr{$target}{"$start..$end"}{$target}{seq} =
            $genome_read_seq;
          $construcTEr{$target}{"$start..$end"}{$te}{start}    = $te_start;
          $construcTEr{$target}{"$start..$end"}{$te}{end}      = $te_end;
          $construcTEr{$target}{"$start..$end"}{$te}{mismatch} = $te_mismatches;
          $construcTEr{$target}{"$start..$end"}{$te}{strand}   = $te_strand;
          $construcTEr{$target}{"$start..$end"}{$te}{id}       = $id;
          $construcTEr{$target}{"$start..$end"}{$te}{seq}      = $te_read_seq;
        }
      }
    } ## end: if ( $te_hit and $genome_aln )
    else {
      delete $seqs{$te}{$id};
    }
  }
}

foreach my $target ( sort keys %construcTEr ) {
  foreach my $range_str (
    sort { ( split /\.\./, $a )[0] <=> ( split /\.\./, $b )[0] }
    keys %{ $construcTEr{$target} }
    )
  {
    foreach my $name ( keys %{ $construcTEr{$target}{$range_str} } ) {
      my $start    = $construcTEr{$target}{$range_str}{$name}{start};
      my $end      = $construcTEr{$target}{$range_str}{$name}{end};
      my $strand   = $construcTEr{$target}{$range_str}{$name}{strand};
      my $seq      = $construcTEr{$target}{$range_str}{$name}{seq};
      my $mismatch = $construcTEr{$target}{$range_str}{$name}{mismatch};
      my $id       = $construcTEr{$target}{$range_str}{$name}{id};

      #if ($strand eq '-'){
      #  ( $seq = reverse $seq ) =~ tr/AaGgTtCcNn/TtCcAaGgNn/;
      #  $id = $id.".revcom";
      #}
      ## if record is a TE hit - get range info

      my $range_str = '';
      if ( exists $ranges{$name} ) {
        my $read_range = range_create( $start, $end );
        foreach my $range_name ( sort keys %{ $ranges{$name} } ) {
          my $te_range = $ranges{$name}{$range_name};
          my $overlap = range_overlap( $te_range, $read_range );
          if ( $overlap >= 5 ) {
            $range_str = $range_str . ",overlap_with_$range_name=$overlap";
          }
        }
        $range_str =~ s/^,//;
      }

      print
">$id $name:$start..$end ($strand) mismatches=$mismatch $range_str\n$seq\n";
    }
  }
}

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

sub range_create {
## range needs to be s<=e
## range is in 1 base notation
## takes 2 numbers and returns an anonymous array
  my @range = ( $_[0], $_[1] );
  if ( $range[0] !~ /^\d+$/ or $range[1] !~ /^\d+$/ ) {
    die "Ranges provided in the TE fasta must be numbers but \n";
  }
  elsif ( $range[0] == 0 or $range[1] == 0 ) {
    die "Ranges are in 1 base notation not 0, numbers must be > 0\n";
  }
  @range = sort { $a <=> $b } @range;
  return \@range;
}

sub range_check {
  my $range = shift;
  if ( ref $range !~ /ARRAY/ ) {
    die "range functions need to be given an array ref\n";
  }
  if ( scalar @$range != 2 ) {
    die "range functions need to have 2 values\n";
  }

  return $range;

}

sub range_overlap {
  ## returns the size of the overlap
  my ( $r1, $r2 ) = @_;
  $r1 = range_check($r1);
  $r2 = range_check($r2);

  my ( $s,  $e )  = @$r1;
  my ( $s2, $e2 ) = @$r2;
  ## r     |-----------|
  ## r2 |-->
  ## r2    |
  if ( $s2 <= $s and $e2 >= $s ) {
    if ( $e2 <= $e ) {
      return ( $e2 - $s + 1 );
    }
    elsif ( $e2 > $e ) {
      return ( $e - $s + 1 );
    }
    else {
      die "range: error1\n";
    }
  }
  ## r  |---------------|
  ## r2 |--->
  ## r2    |--->
  elsif ( $s2 >= $s and $s2 <= $e and $e2 >= $s ) {
    if ( $e2 <= $e ) {
      return ( $e2 - $s + 1 );
    }
    elsif ( $e2 > $e ) {
      return ( $e - $s + 1 );
    }
    else {
      die "range: error2\n";
    }

  }
  else {
    return 0;
  }
}
