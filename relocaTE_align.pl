#!/usr/bin/perl -w
use strict;

#line from relocaTE.pl that calls this script
#$scripts/relocaTE_align.pl $scripts $param_path $genome_file $outregex $TE $exper
my $scripts     = shift;    #full path to scripts directory
my $path        = shift;    #current/top/TE
my $genome_file = shift;
my $regex_file  = shift;
my $TE          = shift;
my $exper       = shift;

##get the regelar expression patterns for mates and for the TE
##when passed on the command line as an argument, even in single
##quotes I lose special regex characters
open INREGEX, "$regex_file" or die "$!";
my $mate_file_1;
my $mate_file_2;
my $mate_file_unpaired;
my $TSD;

while ( my $line = <INREGEX> ) {
  chomp $line;
  ( $mate_file_1, $mate_file_2, $mate_file_unpaired, $TSD ) = split /\t/, $line;
}

my %flanking_fq;
##mate finder
my @files_1;
my @files_2;
my @files_unpaired;
my @flanking_files = <$path/flanking_seq/*flankingReads.fq>;
foreach my $file (@flanking_files) {
  next if -z $file;    ##next file if size is 0
  if ( $file =~ /$mate_file_unpaired/ ) {
    push( @files_unpaired, $file );
  }
  elsif ( $file =~ /$mate_file_1/ ) {
    push( @files_1, $file );
  }
  elsif ( $file =~ /$mate_file_2/ ) {
    push( @files_2, $file );
  }
}
if ( @files_1 and @files_2 ) {
  for ( my $i = 0 ; $i < @files_1 ; $i++ ) {
    my $file_1 = $files_1[$i];
    $file_1 =~ s/$mate_file_1//;
    for ( my $j = 0 ; $j < @files_2 ; $j++ ) {
      my $file_2 = $files_2[$j];
      $file_2 =~ s/$mate_file_2//;
      if ( $file_1 eq $file_2 ) {
        $flanking_fq{$file_1}{1} = $files_1[$i];
        $flanking_fq{$file_1}{2} = $files_2[$j];
        if (@files_unpaired) {
          for ( my $k = 0 ; $k < @files_unpaired ; $k++ ) {
            my $file_unpaired = $files_unpaired[$k];
            if ( $file_1 eq $file_unpaired ) {
              $flanking_fq{$file_1}{unpaired} = $files_unpaired[$k];
              last;
            }
          }
        }
        #if $file_1 eq $file_2 & are finished with unpaired go back to $i loop
        last;
      }
    }
  }
}
else {    ##if only unmatched files are provided
  my @files_singles = <$path/flanking_seq/*flankingReads.fq>;
  foreach my $file ( sort @files_singles ) {
    $flanking_fq{$file}{unpaired} = $file;
  }
}
##get directory path and file name in separate variables
my @genome_dir = split '/', $genome_file;
my $genome_fa  = pop @genome_dir;
my $genome_dir = join '/', @genome_dir;
$genome_fa =~ /(.+)\.fa$/;
my $target     = $1;
my $target_dir = "$path/$target";
##make new directories for file output
`mkdir -p $target_dir`;
`mkdir -p $path/$target/bowtie_aln`;
my $te_dir_path = $path;
$path = $target_dir;
my @bowtie_out_files;

##align each newly created flanking fq files to the genome
##align all files indvidually
##then align again as mates
##will remove any redundant alignments
foreach my $key ( sort keys %flanking_fq ) {
  foreach my $type ( sort keys %{ $flanking_fq{$key} } ) {
    my $flanking_fq = $flanking_fq{$key}{$type};

    #remove and save filename part of path
    my @fq_path = split '/', $flanking_fq;
    my $fq_name = pop @fq_path;
    $fq_name =~ s/\.fq$//;
`bowtie --best -q $genome_file.bowtie_build_index $flanking_fq  1> $path/bowtie_aln/$target.$fq_name.bowtie.single.out 2>> $path/$target.stderr`;
    push @bowtie_out_files,
      "$path/bowtie_aln/$target.$fq_name.bowtie.single.out";
  }    #end of foreach my $type ( sort keys %{ $flanking_fq{$key} } )
  if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2} ) {
    my $flanking_fq_1 = $flanking_fq{$key}{1};
    my $flanking_fq_2 = $flanking_fq{$key}{2};
    my @fq_path       = split '/', $flanking_fq_1;
    my $fq_name       = pop @fq_path;
    $fq_name =~ s/\.fq$//;
    if ( -s $flanking_fq_1 and -s $flanking_fq_2 ) {
      #clean reads if both flanking.fq are non-zero file size
      #clean means make sure the mate1&2 files are in the same order
      #and that any unmated reads are in the unpaired file
`$scripts/clean_pairs_memory.pl -1 $flanking_fq_1 -2 $flanking_fq_2 1> $te_dir_path/flanking_seq/$fq_name.unPaired.fq 2>> $path/$target.stderr`;
    }    #end of if ( -s $flanking_fq_1 and -s $flanking_fq_2 )
    if (  -s "$flanking_fq_1.matched"
      and -s "$flanking_fq_2.matched" )
    {
`bowtie --best -q $genome_file.bowtie_build_index -1 $flanking_fq_1.matched -2 $flanking_fq_2.matched 1> $path/bowtie_aln/$target.$fq_name.bowtie.mates.out 2>> $path/$target.stderr`;
      push @bowtie_out_files,
        "$path/bowtie_aln/$target.$fq_name.bowtie.mates.out";
`bowtie --best -q $genome_file.bowtie_build_index $te_dir_path/flanking_seq/$fq_name.unPaired.fq 1> $path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out 2>> $path/$target.stderr`;
      push @bowtie_out_files,
        "$path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out";
    }  # end of if(-s "$flanking_fq_1.matched" and -s "$flanking_fq_2.matched" )
  }  #end of if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2})
}    #end of foreach $key
my @files2merge;
foreach my $bowtie_out (@bowtie_out_files) {
  if ( -s $bowtie_out ) {
    #covert bowtie output to sam
    #if using the bowtie option to produce sam files i am unable to get rid of duplicate alignments due to flags in col2
    `bowtie2sam.pl $bowtie_out 1> $bowtie_out.sam 2>> $path/$target.stderr`;

    #convert sam to bam
`samtools import $genome_dir/$genome_fa.fai $bowtie_out.sam $bowtie_out.bam 2>> $path/$target.stderr`;

    #sort bam
    `samtools sort $bowtie_out.bam $bowtie_out.sorted 2>> $path/$target.stderr`;

    #index bam
    `samtools index $bowtie_out.sorted.bam 2>> $path/$target.stderr`;
    push @files2merge, "$bowtie_out.sorted.bam" if -s "$bowtie_out.sorted.bam";
  }    #end of if ( -s $bowtie_out )
}    #end of foreach my $bowtie_out (@bowtie_out_files)
my $files2merge;
##if there are too many bam files to merge samtools can't handle it
##so make smaller batches of merging,then merge the tmp bam files
my $bam_file_count = scalar @files2merge;
if ( @files2merge > 50 ) {
  my $filecount = scalar @files2merge;
  my @big_files_2_merge;
  for ( my $i = 0 ; $i < $filecount ; $i = $i + 50 ) {
    $files2merge = join " ", ( splice( @files2merge, 0, 50 ) );
    if (scalar @files2merge > 1){
`samtools merge -f $path/bowtie_aln/$target.$TE.merged.bam.$i.temp $files2merge 2>> $path/$target.stderr`;
      push @big_files_2_merge, "$path/bowtie_aln/$target.$TE.merged.bam.$i.temp";
    }else {
      ##if there is only 1 file after processing the 50, then just push onto @big_files_2_merge
      push @big_files_2_merge, $files2merge;
    }
  }
  $files2merge = join ' ', @big_files_2_merge if @big_files_2_merge;
}
else {
  $files2merge = join ' ', @files2merge if @files2merge;
}
#merge paired and unPaired bam
if (defined $files2merge){
   my $bam_file;
   if ($bam_file_count > 1){
    `samtools merge -f $path/bowtie_aln/$target.$TE.merged.bam $files2merge 2>> $path/$target.stderr`;
    `samtools sort $path/bowtie_aln/$target.$TE.merged.bam $path/bowtie_aln/$target.$TE.merged.sorted 2>> $path/$target.stderr`;
  `  samtools index $path/bowtie_aln/$target.$TE.merged.sorted.bam 2>> $path/$target.stderr`;
     $bam_file = "$path/bowtie_aln/$target.$TE.merged.sorted.bam";
   }else { ##only one file. can't merge 1 file.
     $bam_file = $files2merge;
   }
  #identify mping insertion sites
  `$scripts/relocaTE_insertionFinder.pl $bam_file $target $genome_file $TE $TSD $exper`;
}
