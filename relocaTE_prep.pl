#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Cwd;
use strict;

##change $scripts to location of find_TE scripts
my $scripts = '~/bin';
my $current_dir = shift;
my $top_dir = shift;
my $genomeFasta = shift;
my $te_fasta = shift;

##split genome file into individual fasta files
my @genome_fastas;
if ($mapping){
  my $genome_path = File::Spec->rel2abs($genomeFasta);
  open( INFASTA, "$genome_path" ) || die "$!\n";
  my $i = 0;
  while ( my $line = <INFASTA> ) {
    if ( $line =~ /^>(\S+)/ ) {
      my $id = $1;
      $id = s/\W/_/g;
      if ( $i > 0 ) {
        close(OUTFASTA);
        $i = 0;
      }
      my @genome_dir = split '/'  , $genome_path;
      pop @genome_dir;
      $genome_dir = join '/' , @genomeDir;
      my $new_file = "$genome_dir/$id.fa";
      push @genome_fastas, $new_file;
      open( OUTFASTA, ">$new_file" ) or die "$!\n";
      print OUTFASTA $line;
      $i++;
    }
    elsif ($line !~ /^>/) {  ##should be sequence
      print OUTFASTA $line;
    }else {
      die "Your genome fasta file is in a unexpected format. 
I was expecting a line of seqeunce but found something else:
$line\n";
    }
  }
  close(INFASTA);
  close(OUTFASTA);
  foreach my $genome_file (@genome_fastas){
    if ( !-e "$genome_file.bowtie_build_index.1.ebwt" and $mapping ) {
      `bowtie-build -f $genome_file $genome_file.bowtie_build_index`;
    }
    #create an index of genome fasta
    `samtools faidx $genome_file`;
  }
}##end if($mapping)

my $te_path = File::Spec->rel2abs($te_fasta);
#convert fq files to fa for blat
my @fq;
my @fa;

foreach my $fq (@fq_files) {
    my $fq_path = File::Spec->rel2abs($fq);
    push @fq, $fq_path;
    my $fa = $fq;
    if ( $fa =~ s/\.(fq|fastq)$/.fa/ ) {
        push @fa, $fa;
        if ( !-e $fa ) {
            open INFQ,  $fq_path or die $!;
            open OUTFA, ">$fa"   or die $!;

            while ( my $header = <INFQ> ) {
                my $seq         = <INFQ>;
                my $qual_header = <INFQ>;
                my $qual        = <INFQ>;

                die "ERROR: expected \'\@\' but saw $header"
                  if substr( $header, 0, 1 ) ne '@';

                print OUTFA ">", substr( $header, 1 );
                print OUTFA $seq;
            }
            close INFQ;
            close OUTFA;
        }
    }
    else {
        print
"$fq does not seem to be a fastq based on the file extension. It should be fq or fastq\n";
        &getHelp();
    }
}

##split TE fasta into single record fastas
my @te_fastas;
my %TSD;

##put in a new directory structure workingDir/date-te-search/tefilename/all-newly-created-files
##create new te fasta file
my $top_dir = $outdir;
open( INFASTA, "$te_fasta" ) || die "$!\n";
my $i = 0;
while ( my $line = <INFASTA> ) {
    if ( $line =~ /^>(\S+)\s+TSD=(\S+)/ ) {
        my $id = $1;
        $TSD{$id} = $2;
        if ( $i > 0 ) {
            close(OUTFASTA);
            $i = 0;
        }
        my $te_file = "$id.fa";
        $te_file =~ s/\|/_/g;
        my $te_dir = "$current_dir/$top_dir/$id";
        push @te_fastas, "$te_dir/$te_file";
        `mkdir -p $te_dir`;
        open( OUTFASTA, ">$te_dir/$te_file" ) or die "$!\n";
        print OUTFASTA $line;
        $i++;
    }elsif ($line =~/^>/ and $line !~ /TSD=/){
        die  "The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME TSD=TSD\nSEQUENCE\n";
    }
    else {  ##should be sequence
        print OUTFASTA $line;
    }
}
close(INFASTA);
close(OUTFASTA);

#foreach TE fasta blat against target chromosome and parse and find insertion sites
foreach my $te_path (@te_fastas) {
    my @path     = split '/', $te_path;
    my $te_fasta = pop @path;
    my $path     = join '/', @path;
    my $TE       = $te_fasta;
    $TE =~ s/\.fa//;

    #blat fa files against te.fa
    my @flanking_fq;
    my @flanking_fq_mates;
    my $fq_file_count = scalar @fq;
    for ( my $i = 0 ; $i < $fq_file_count ; $i++ ) {
        my $fa = $fa[$i];
        my $fq = $fq[$i];

        #remove and save filename part of path
        my @fa_path = split '/', $fa;
        my $fa_name = pop @fa_path;
        $fa_name =~ s/\.fa$//;
        
        ##create a directory for each genome seq
        $path = "$path/$fa_name";
        `mkdir -p $path`;
`blat -minScore=10 -tileSize=7 $te_path $fa $path/$fa_name.te_$TE.blatout`
          if !-e "$path/$fa_name.te_$TE.blatout";

        #my $file_num         = $i + 1;
        my $te_Containing_fq = "$path/$fa_name.te_$TE.ContainingReads.fq";
        if ( -e $te_Containing_fq ) {
            $fq = $te_Containing_fq;
        }
`perl $scripts/get_fq_of_te_trimmed_te_matching_reads.pl $path/$fa_name.te_$TE.blatout $fq $len_cutoff $mismatch_allowance > $path/$fa_name.te_$TE.flankingReads.fq `;
    }

##insert mate finder here
    my %flanking_fq;
    my @files_1;
    my @files_2;
    my @files_unpaired;
    my @flanking_files = <$path/*flankingReads.fq>;
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
    my @files_singles = <$path/*flankingReads.fq>;
    foreach my $file ( sort @files_singles ) {
      $flanking_fq{$file}{unpaired} = $file;
    }
  }
  die "Problem finding *flankingReads.fq" if scalar (keys %flanking_fq) == 0;
}
open OUTFQFILE , ">$fq_file_list" or die $!;
foreach my $base (keys %flanking_fq){
  foreach my $pair ( keys %{$flanking_fq{$base}}){
    my $file = $flanking_fq{$base}{$pair};
    print OUTFQFILE "$base\t$pair\t$file\n";
  }
}

my $fq_file_list = "$current_dir/$top_dir/fq_files.list";
foreach my $te_file (@te_fastas){
 foreach my $genome_file (@genome_fastas){
  `$scripts/relocaTE_prep.pl $current_dir $top_dir $genome_file $te_file $fq_file_list`;
 }
}
