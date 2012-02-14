#!/usr/bin/perl -w
use File::Spec;
use Getopt::Long;
use Cwd;
use strict;

##change $scripts to location of relocaTE scripts
my $scripts = '~/bin/relocaTE_editing';

if ( !defined @ARGV ) {
  &getHelp();
}
my $genomeFasta = 'NONE';
my $te_fasta;
my $target             = 'NONE';
my $len_cutoff         = 10;
my $mismatch_allowance = 0;
my $fq_dir;
my $exper = 'not.given';
## REGEX = _1, followed by (\D*?) optional non-digit characters, followed by fq
my $mate_file_1        = '_1\D*?fq';
my $mate_file_2        = '_2\D*?fq';
my $mate_file_unpaired = '.unPaired\D*?fq';
my $workingdir;
my $outdir   = 'outdir_teSearch';
my $parallel = 1;
GetOptions(
  'p|parallel:i'    => \$parallel,
  'e|exper:s'       => \$exper,
  'w|workingdir:s'  => \$workingdir,
  'o|outdir:s'      => \$outdir,
  'd|fq_dir:s'      => \$fq_dir,
  'g|genomeFasta:s' => \$genomeFasta,
  't|te_fasta:s'    => \$te_fasta,
  'l|len_cutoff:i'  => \$len_cutoff,
  'm|mismatch:f'    => \$mismatch_allowance,
  '1|mate_1_id:s'   => \$mate_file_1,
  '2|mate_2_id:s'   => \$mate_file_2,
  'u|unpaired_id:s' => \$mate_file_unpaired,
  'h|help'          => \&getHelp,
);
my $current_dir;

if ( defined $workingdir and -d $workingdir ) {
  $current_dir = File::Spec->rel2abs($workingdir);
  $current_dir =~ s/\/$//;
}
else {
  $current_dir = cwd();
}
my $mapping = 1;

if ( !defined $genomeFasta ) {
  print "\n\nPlease provide reference genome by using -g Genome fasta path\n";
  &getHelp();
}
elsif ( $genomeFasta eq 'NONE' ) {
  print
"You did not provide a genome fasta, if you proceed only reads containing the TE will be found, no mapping of insertions will be performed\n";
  print "Proceed without mapping?\n";
  my $answer;
  while ( $answer = <STDIN> ) {

    #Exit if it was just spaces (or just an enter)
    last if $answer =~ /^\s*|\n$/;
  }
  if ( $answer =~ /n/i ) {
    &getHelp();
  }
  else {
    $mapping = 0;
  }
}
elsif ( !-e $genomeFasta ) {
  print "$genomeFasta does not exist. Check file name.\n";
  &getHelp();
}
if ( !defined $te_fasta ) {
  print
"\n\nPlease provide fasta file containing transposable elements by using -t TE fasta path\n";
  &getHelp();
}
elsif ( !-e $te_fasta ) {
  print "$te_fasta does not exist. Check file name.\n";
  &getHelp();
}
else {
  my $first_line = `head -n1 $te_fasta`;
  if ( $first_line !~ /^>\S+\s+TSD=\S+/ ) {
    die
"The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME TSD=TSD\nSEQUENCE\n";
  }
}
my @fq_files;
my %fq_files;
if ( !defined $fq_dir ) {
  print "\n\nPlease provide a directory of paired fastq files\n";
  &getHelp();
}
elsif ( !-d $fq_dir ) {
  print
"\n\nCheck the spelling or location of $fq_dir, Please provide a directory of paired fastq files\n";
  &getHelp();
}
else {
  my $fq_path = File::Spec->rel2abs($fq_dir);
  @fq_files = <$fq_path/*fq>;
  my @fastq_files = <$fq_path/*fastq>;
  push @fq_files, @fastq_files;
  if ( scalar @fq_files == 0 ) {
    print "Must provide at least 1 short read file\n";
    &getHelp();
  }
}

sub getHelp {
  print ' 
usage:
./relocaTE.pl [-t TE_fasta_file][-g chromosome_genome_fasta][-d dir_of_fq][-e short_sample_name][-h] 

options:

**required:
-g STR          genome fasta file path. If not provided will only align reads to TE and remove TE seq from short reads. [no default]
-t STR          fasta containing nucleotide sequences of transposable elements with TSD=xxx in the desc. [no default]
-d STR          directory of paired and unpaired fastq files (paired _1.fq & _2.fq) (.fq or .fastq is acceptable)  [no default]

**recommended:
-e STR          Short sample name, will be used in the output files to create IDs for the insert (ex. A123) [not.given]
-o STR          name for directory to contain output directories and files, will be created for the run (ex. 04222012_A123) [outdir_teSearch]

**optional:
-p INT          run each genome sequence separetly, parallel. The alternative (0) would be to run one after the other (int, 0=false or 1=true) [1] 
-w STR          base working directory, needs to exist, will not create, full path [cwd] 
-l INT          len cutoff for the te trimmed reads to be aligned [10] 
-m FRACTION     mismatch allowance for alignment to TE (int, ex 0.1) [0] 
-1 STR		regular expression to identify mate 1 paired files [_1\D*?fq]
-2 STR          regular expression to identify mate 2 paired files [_2\D*?fq]
-u STR          regular expression to identify unpaied files [.unPaired\D*?fq] 
-h              this message

SAMPLE TE FASTA
>mping TSD=TTA
GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATG
ATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTT
TCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGT
CCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGGATTTTGGGTACAAATGATCCCAGCAA
CTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGT
TTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATT
GTGACTGGCC

';

  exit 1;
}
$outdir =~ s/\/$//;
my $te_path = File::Spec->rel2abs($te_fasta);
my $top_dir = $outdir;

##split genome file into individual fasta files
my @genome_fastas;
if ($mapping) {
  my $genome_path = File::Spec->rel2abs($genomeFasta);
  open( INFASTA, "$genome_path" ) || die "$!\n";
  my $i = 0;
  my $exists = 0;
  while ( my $line = <INFASTA> ) {
    if ($exists){
      next unless $line =~ /^>(\S+)/;
    }
    if ( $line =~ /^>(\S+)/ ) {
      my $id = $1;
      if ( $id =~ /\|/ ) {
        $id   =~ s/\|/_/g;
        $line =~ s/\|/_/g;
      }
      if ( $i > 0 ) {
        close(OUTFASTA);
        $i = 0;
      }
      my @genome_dir = split '/', $genome_path;
      pop @genome_dir;
      my $genome_dir = join '/', @genome_dir;
      my $new_file = "$genome_dir/$id.fa";
      push @genome_fastas, $new_file;
      if (-e $new_file){
        $exists = 1;
        next;
      }else{
        $exists = 0;
      }
      open( OUTFASTA, ">$new_file" ) or die "$!\n";
      print OUTFASTA $line;
      $i++;
    }
    elsif ( $line !~ /^>/ ) {    ##should be sequence
      print OUTFASTA $line;
    }
    else {
      die "Your genome fasta file is in a unexpected format. 
I was expecting a line of seqeunce but found something else:
$line\n";
    }
  }
  close(INFASTA);
  close(OUTFASTA);
  ##format genome sequences for bowtie and samtools
  foreach my $genome_file (@genome_fastas) {
    my @cmds;

    #create bowtie index
    if ( !-e "$genome_file.bowtie_build_index.1.ebwt" ) {
      push @cmds,
        "bowtie-build -f $genome_file $genome_file.bowtie_build_index";
    }

    #create an index of genome fasta
    if ( !-e "$genome_file.fai" ) {
      push @cmds, "samtools faidx $genome_file";
    }
    $genome_file =~ /.+\/(.+)\.fa$/;
    my $ref = $1;
    if ( $parallel and @cmds ) {
      my $shell_dir = "$current_dir/$top_dir/shellscripts/step_1";
      `mkdir -p $shell_dir`;
      open OUTSH, ">$shell_dir/$ref.formatGenome.sh";
      my $cmd = join "\n", @cmds;
      print OUTSH "$cmd\n";
      close OUTSH;
    }
    elsif ( $parallel and !@cmds ) {
      my $step2_file =
"$current_dir/$top_dir/shellscripts/step_1_not_needed_genomefasta_already_formatted";
      my $shell_dir = "$current_dir/$top_dir/shellscripts";
      `mkdir -p $shell_dir`;
      `touch $step2_file` if $parallel;
    }
    elsif (@cmds) {
      ##run them serierally now
      foreach my $cmd (@cmds) {
        `$cmd`;
      }
    }
  }
}    ##end if($mapping)

my @fq;
my @fa;

#convert fq files to fa for blat
foreach my $fq (@fq_files) {
  my $fq_path = File::Spec->rel2abs($fq);
  push @fq, $fq_path;
  my $fa = $fq;
  if ( $fa =~ s/\.(fq|fastq)$/.fa/ ) {
    push @fa, $fa;
    if ( !-e $fa ) {
      my $cmd = "$scripts/relocaTE_fq2fa.pl $fq_path $fa";
      if ($parallel) {
        my @fq_path   = split '/', $fq_path;
        my $fq_name   = pop @fq_path;
        my $shell_dir = "$current_dir/$top_dir/shellscripts/step_2";
        `mkdir -p $shell_dir`;
        my $outsh = ">$shell_dir/$fq_name" . "2fa.sh";
        open OUTSH, ">$outsh";
        print OUTSH "$cmd\n";
      }
      else {
        `$cmd`;
      }
    }
    else {
      my $shell_dir = "$current_dir/$top_dir/shellscripts";
      `mkdir -p $shell_dir`;
      my $step2_file =
"$current_dir/$top_dir/shellscripts/step_2_not_needed_fq_already_converted_2_fa";
      `touch $step2_file` if $parallel;
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

##split the TE fasta of many seqs into individual files
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
    ##create new dir for files: workingDir/outdir/TE/
    my $te_dir = "$current_dir/$top_dir/$id";
    push @te_fastas, "$te_dir/$te_file";
    `mkdir -p $te_dir`;
    open( OUTFASTA, ">$te_dir/$te_file" ) or die "$!\n";
    print OUTFASTA $line;
    $i++;
  }
  elsif ( $line =~ /^>/ and $line !~ /TSD=/ ) {
    die
"The TE_fasta:$te_fasta does not have the proper format:\n>TE_NAME TSD=TSD\nSEQUENCE\n";
  }
  else {    ##should be sequence
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
  `mkdir -p $path/blat_output`;
  `mkdir -p $path/flanking_seq`;
  `mkdir -p $path/te_containing_fq`;
  `mkdir -p $path/te_only_read_portions_fa`;

  #blat fa files against te.fa
  my @flanking_fq;
  my $fq_file_count = scalar @fq;
  for ( my $i = 0 ; $i < $fq_file_count ; $i++ ) {
    my $fa = $fa[$i];
    my $fq = $fq[$i];

    #remove and save filename part of path
    my @fa_path = split '/', $fa;
    my $fa_name = pop @fa_path;
    $fa_name =~ s/\.fa$//;
    if ($parallel) {
      my $shell_dir = "$current_dir/$top_dir/shellscripts/step_3/$TE";
      `mkdir -p $shell_dir`;
      open OUTSH, ">$shell_dir/$TE.$fa_name.blat.sh";
    }
    if ( !-e "$path/blat_output/$fa_name.te_$TE.blatout" ) {
      my $cmd =
"blat -minScore=10 -tileSize=7 $te_path $fa $path/blat_output/$fa_name.te_$TE.blatout";
      print OUTSH "$cmd\n" if $parallel;
      `$cmd` if !$parallel;
    }
    my $te_Containing_fq =
      "$path/te_containing_fq/$fa_name.te_$TE.ContainingReads.fq";
    if ( -e $te_Containing_fq ) {
      $fq = $te_Containing_fq;
    }
    my $cmd =
"perl $scripts/relocaTE_trim.pl $path/blat_output/$fa_name.te_$TE.blatout $fq $len_cutoff $mismatch_allowance > $path/flanking_seq/$fa_name.te_$TE.flankingReads.fq";
    if ($parallel) {
      print OUTSH "$cmd\n";
      close OUTSH;
    }
    else {
      `$cmd` if !$parallel;
    }
  }
  ##if a genome file was provided, align seqs to genome
  ##if no genome file was provided, will only blat and trim reads of te seq
  if ($mapping) {
    foreach my $genome_file (@genome_fastas) {
      my $param_path = "$current_dir/$top_dir/$TE";
      my $outregex   = "$param_path/regex.txt";
      open OUTREGEX, ">$outregex" or die $!;
      print OUTREGEX
        "$mate_file_1\t$mate_file_2\t$mate_file_unpaired\t$TSD{$TE}";
      my $cmd =
"$scripts/relocaTE_align.pl $scripts $param_path $genome_file $outregex $TE $exper";
      if ( !$parallel ) {
        `$cmd`;
      }
      else {
        my $shell_dir = "$current_dir/$top_dir/shellscripts/step_4/$TE";
        $genome_file =~ /.+\/(.+)\.fa$/;
        my $ref = $1;
        `mkdir -p $shell_dir`;
        open OUTSH, ">$shell_dir/$TE.$ref.findSites.sh";
        print OUTSH "$cmd\n";
        close OUTSH;
      }
    }
  }
}
