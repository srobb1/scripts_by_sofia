#!/usr/bin/perl -w
use strict;
use Data::Dumper;

#$scripts/relocaTE_process.pl $scripts $path $genome_file $mate_file_1 $mate_file_2 $mate_file_unpaired $TE $TSD{$TE} $exper
my $scripts             = shift;    #full path to scripts directory
my $path                = shift;    #current/top/TE
my $genome_file         = shift;
my $regex_file          = shift;
my $TE                  = shift;
my $exper               = shift;

open INREGEX , "$regex_file" or die "$!";
my $mate_file_1;   
my $mate_file_2;    
my $mate_file_unpaired;
my $TSD;
while (my $line = <INREGEX>){
  chomp $line;
  ($mate_file_1,$mate_file_2,$mate_file_unpaired,$TSD) = split /\t/ , $line ;
}

my %flanking_fq;
##insert mate finder here
    my @files_1;
    my @files_2;
    my @files_unpaired;
    my @flanking_files = <$path/flanking_seq/*flankingReads.fq>;
    foreach my $file  (@flanking_files){
      next if -z $file; ##next file if size is 0
      if ($file =~ /$mate_file_unpaired/){
        push (@files_unpaired , $file);
      }elsif ($file =~ /$mate_file_1/){
        push (@files_1, $file);
      }elsif ($file =~ /$mate_file_2/){
        push (@files_2, $file);
      }
    }
#print "1:$mate_file_1 ", scalar @files_1 , " 2:$mate_file_2 " , scalar @files_2 , " un:$mate_file_unpaired " , scalar @files_unpaired ,"\n";
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
                        $flanking_fq{$file_1}{unpaired} =
                          $files_unpaired[$k];
                        last;
                    }
                }
            }
            #if $file_1 eq $file_2 and we are finished with unpaired go back to $i loop
            last;
          }
       }
    }
  }
  else {              ##if only unmatched files are provided
    my @files_singles = <$path/flanking_seq/*flankingReads.fq>;
    foreach my $file ( sort @files_singles ) {
      $flanking_fq{$file}{unpaired} = $file;
    }
  }
#print Dumper \%flanking_fq;
my @genome_dir = split '/' , $genome_file;
my $genome_fa        = pop @genome_dir;
my $genome_dir = join '/' , @genome_dir;
$genome_fa =~ /(.+)\.fa$/;
my $target     = $1;
my $target_dir = "$path/$target";
`mkdir -p $target_dir`;
`mkdir -p $path/$target/bowtie_aln`;
my $te_dir_path = $path;
$path = $target_dir;
my @bowtie_out_files;

foreach my $key ( sort keys %flanking_fq ) {
    foreach my $type ( sort keys %{ $flanking_fq{$key} } ) {
        my $flanking_fq = $flanking_fq{$key}{$type};

        #remove and save filename part of path
        my @fq_path = split '/', $flanking_fq;
        my $fq_name = pop @fq_path;
        $fq_name =~ s/\.fq$//;
`bowtie --best -q $genome_file.bowtie_build_index $flanking_fq  1> $path/bowtie_aln/$target.$fq_name.bowtie.single.out 2>> $path/$target.stderr`;
        push @bowtie_out_files, "$path/bowtie_aln/$target.$fq_name.bowtie.single.out";
    }    #end of foreach my $type ( sort keys %{ $flanking_fq{$key} } )
    if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2} ) {
        my $flanking_fq_1 = $flanking_fq{$key}{1};
        my $flanking_fq_2 = $flanking_fq{$key}{2};
        my @fq_path       = split '/', $flanking_fq_1;
        my $fq_name       = pop @fq_path;
        $fq_name =~ s/\.fq$//;
        if ( -s $flanking_fq_1 and -s $flanking_fq_2 ) {

            #clean reads if both flanking.fq are non-zero file size
`$scripts/clean_pairs_memory.pl -1 $flanking_fq_1 -2 $flanking_fq_2 1> $te_dir_path/flanking_seq/$fq_name.unPaired.fq 2>> $path/$target.stderr`;
        }    #end of if ( -s $flanking_fq_1 and -s $flanking_fq_2 )
        if (    -s "$flanking_fq_1.matched"
            and -s "$flanking_fq_2.matched" )
        {
`bowtie --best -q $genome_file.bowtie_build_index -1 $flanking_fq_1.matched -2 $flanking_fq_2.matched 1> $path/bowtie_aln/$target.$fq_name.bowtie.mates.out 2>> $path/$target.stderr`;
            push @bowtie_out_files, "$path/bowtie_aln/$target.$fq_name.bowtie.mates.out";
`bowtie --best -q $genome_file.bowtie_build_index $te_dir_path/flanking_seq/$fq_name.unPaired.fq 1> $path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out 2>> $path/$target.stderr`;
            push @bowtie_out_files,
              "$path/bowtie_aln/$target.$fq_name.bowtie.unPaired.out";
        } # end of if(-s "$flanking_fq_1.matched" and -s "$flanking_fq_2.matched" )
    } #end of if ( exists $flanking_fq{$key}{1} and exists $flanking_fq{$key}{2})
}    #end of foreach $key
my @files2merge;
foreach my $bowtie_out (@bowtie_out_files) {
    if ( -s $bowtie_out ) {
        #covert bowtie output to sam
        `bowtie2sam.pl $bowtie_out 1> $bowtie_out.sam 2>> $path/$target.stderr`;

        #convert sam to bam
`samtools import $genome_dir/$genome_fa.fai $bowtie_out.sam $bowtie_out.bam 2>> $path/$target.stderr`;
#print "samtools import $genome_dir/$genome_fa.fa.fai $bowtie_out.sam $bowtie_out.bam 2>> $path/$target.stderr\n";

        #sort bam
`samtools sort $bowtie_out.bam $bowtie_out.sorted 2>> $path/$target.stderr`;
#print "samtools sort $bowtie_out.bam $bowtie_out.sorted 2>> $path/$target.stderr\n";

        #index bam
        `samtools index $bowtie_out.sorted.bam 2>> $path/$target.stderr`;
        push @files2merge, "$bowtie_out.sorted.bam 2>> $path/$target.stderr";
    }    #end of if ( -s $bowtie_out )
}    #end of foreach my $bowtie_out (@bowtie_out_files)
my $files2merge;
if ( @files2merge > 50 ) {
    my $filecount = scalar @files2merge;
    my @big_files_2_merge;
    for ( my $i = 0 ; $i < $filecount ; $i = $i + 50 ) {
        $files2merge = join " ", ( splice( @files2merge, 0, 50 ) );
`samtools merge -f $path/bowtie_aln/$target.$TE.merged.bam.$i.temp $files2merge 2>> $path/$target.stderr`;
        push @big_files_2_merge, "$path/bowtie_aln/$target.$TE.merged.bam.$i.temp";
    }

    #$files2merge = join ', ' , @big_files_2_merge;
    $files2merge = "$path/bowtie_aln/$target.$TE.merged.bam.*.temp";
}
else {
    $files2merge = join ', ', @files2merge;
}

#merge paired and unPaired bam
`samtools merge -f $path/bowtie_aln/$target.$TE.merged.bam $files2merge 2>> $path/$target.stderr`;
`samtools sort $path/bowtie_aln/$target.$TE.merged.bam $path/bowtie_aln/$target.$TE.merged.sorted 2>> $path/$target.stderr`;
`samtools index $path/bowtie_aln/$target.$TE.merged.sorted.bam 2>> $path/$target.stderr`;

#identify mping insertion sites
`$scripts/relocaTE_insertionFinder.pl $path/bowtie_aln/$target.$TE.merged.sorted.bam $target $genome_file $TE $TSD $exper`;
