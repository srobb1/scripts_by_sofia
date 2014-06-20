#!/usr/bin/perl -w
use strict;
use Bio::DB::Fasta;
use Bio::SeqIO;

## need refFasta, TEFasta, RelocaTEouput, padding
my $dbfile     = shift;
my $TEfile     = shift;
my $range_file = shift;
my $padding    = shift;
my $usr_TE; #to be determined from range_file = 'mping';#shift; ## ping or mping
open GFF, ">$range_file.gff" or die "can't open $range_file.gff for writing\n";
if ( !defined $padding ) {
  $padding = 0;
}

if ( !defined $dbfile or !defined $range_file ) {
  die "Please provide fasta file and RelocaTE output\n
example:
./create_any_complete_site_with_flanking.pl RefFastaFile TEFastaFile location_File padding
"
}

my $db_obj = Bio::DB::Fasta->new($dbfile);

#my $db_TE_obj = Bio::DB::Fasta->new($TEfile);

my %TEs;
my $seqIO_obj = Bio::SeqIO->new( -format => 'fasta', -file => $TEfile );
while ( my $seq_obj = $seqIO_obj->next_seq ) {
  my $TE  = $seq_obj->id;
  my $seq = $seq_obj->seq;
  $TEs{$TE} = $seq;
}
my %inserts;
open IN, $range_file or die "Can't open range_file\n";
while ( my $range = <IN> ) {
  chomp $range;
  my @range = split /\t/, $range;
  next if $range =~ /^TE.+insertion_site/;

#TE      TSD     Exper   insertion_site  strand  left_flanking_read_count        right_flanking_read_count       left_flanking_seq       right_flanking_seq
#Dasheng ATCCA   A119_2  Chr1:861881..861885     -       32      28      CCACAGAAATTACAACAGCCAAACAACACATGTGGCCAATTTTCACACATAACGAGATAACATCCCAATCACAACTTTATGATAAACGAAAGCTTCCATC    CAGTTCACGGCAAAGAGATTACTCATCATAAATAAATAGTGTACAATAAAATATATCACTGTTTACAAATATTATTTCTAGCTTCAGCATTTTACCACCT
  my ( $ref, $start, $end ) = $range[3] =~ /^(\S+):(\d+)\.\.(\d+)/;
  $inserts{$ref}{$start}{$end}{strains}{ $range[2] }++;
  $inserts{$ref}{$start}{$end}{strand}{ $range[4] }++;

  #$inserts{$ref}{$start}{$end}{TSD} = $range[1];
  $inserts{$ref}{$start}{$end}{TE}{ $range[0] } = $range[1];    ##TSD info
}
foreach my $ref ( keys %inserts ) {
  foreach my $start ( keys %{ $inserts{$ref} } ) {
    foreach my $end ( keys %{ $inserts{$ref}{$start} } ) {
      my @TEs;
      foreach my $TE_name ( keys %{ $inserts{$ref}{$start}{$end}{TE} } ) {
        my $TSD     = $inserts{$ref}{$start}{$end}{TE}{$TE_name};
        my @strains = sort keys %{ $inserts{$ref}{$start}{$end}{strains} };
        push @TEs, $TE_name;
        my $desc = join( "/", @strains );
        my $plus =
          defined $inserts{$ref}{$start}{$end}{strand}{'+'}
          ? $inserts{$ref}{$start}{$end}{strand}{'+'}
          : 0;
        my $minus =
          defined $inserts{$ref}{$start}{$end}{strand}{'-'}
          ? $inserts{$ref}{$start}{$end}{strand}{'-'}
          : 0;
        my $strand_ref =
          defined $inserts{$ref}{$start}{$end}{strand}{'.'}
          ? $inserts{$ref}{$start}{$end}{strand}{'.'}
          : 0;
        my $strand = 'fwd';

        if ( $minus > $plus ) {
          $strand = 'rc';
        } elsif ( $strand_ref > $minus and $strand_ref > $plus ) {
          $strand = 'ref';
        }
        my $name = $TE_name . "." . join( '_', $ref, $start, $end );
        my $seq_obj = $db_obj->get_Seq_by_id($ref);
        if ($seq_obj) {
          my $seq = $seq_obj->seq;
          my %coords;
          my ( $padded_start, $padded_end, $pre_TE, $post_TE, $justTE );
          my $rev_TE;
          my $revcomp_TE;
          my $revcomp_rev_TE;
          if ( defined $seq ) {
            #Chr1:2640500..2640502 Ping1 shared in H/E/19/23 and is in a stowawy
            if ($TE_name eq 'ping' and $ref eq 'Chr1' and $start == 2640500 and $end == 2640502){
              print_special_ping($start,$end,$seq,$name);
              next;
            }
            my $TE = substr $seq, $start - 1, ( $end - $start ) + 1;

            #if this is a reference mPing, add the TSD inaddition to flanking
            if ( length $TE > length $TSD or $TSD eq 'TSD') {
              $justTE= $TE;
              $padded_start = $start - $padding;
              $padded_end   = $end + $padding;
              $pre_TE       = substr $seq, $padded_start - 1, $padding;
              $post_TE      = substr $seq, $end, ($padding);
            }

            ## if TE length is 3, then we need to add the TE and dup the TSD
            elsif ( length $TE == length $TSD ) {
              $TSD = $TE; ## reset the $TSD to be exactly what is in the genomeSeq, just incase
              $justTE = $TEs{$TE_name};
              my $rev_justTE     = reverse $justTE;
              my $revcomp_justTE = $rev_justTE;
              $revcomp_justTE =~ tr/ATGC/TACG/;

              $TE           = $justTE . $TSD;
              $revcomp_TE   = $revcomp_justTE . $TSD;
              $padded_start = $end - $padding;
              $padded_end   = $end + $padding - ( length $TSD ) - 1;
              $pre_TE       = substr $seq, $end - $padding, $padding;
              $post_TE      = substr $seq, $end, $padding - ( length $TSD );
            }
            my $subseq = $pre_TE . $TE . $post_TE;
            ## coords are relative to subseq, always starts with 1 and ends at length of substr
            $coords{flank_1_start} = 1;
            $coords{flank_1_end}   = $padding;
            $coords{TE_start}      = $padding + 1;
            $coords{TE_end}        = $padding + ( length $justTE );
            $coords{TSD_2_start}   = $TSD ne 'TSD' ? $coords{TE_end} + 1 : 'ref';
            $coords{TSD_2_end}     = $TSD ne 'TSD' ? $coords{TE_end} + ( length $TSD ) : 'ref';
            $coords{TSD_1_start}   = $TSD ne 'TSD' ? $coords{TE_start} - ( length $TSD ) : 'ref';
            $coords{TSD_1_end}     = $TSD ne 'TSD' ? $coords{TE_start} - 1 : 'ref' ;
            $coords{flank_2_start} = $coords{TE_end} + 1;
            $coords{flank_2_end}   = $padding * 2 + ( length $justTE );
            my $coords =
              "FLANK1:$coords{flank_1_start}..$coords{flank_1_end} TSD1:$coords{TSD_1_start}..$coords{TSD_1_end} TE:$coords{TE_start}..$coords{TE_end} TSD2:$coords{TSD_2_start}..$coords{TSD_2_end} FLANK2:$coords{flank_2_start}..$coords{flank_2_end}";
            my $revcomp_subseq;

            if ( defined $revcomp_TE ) {
              $revcomp_subseq = $pre_TE . $revcomp_TE . $post_TE;
            }
            my $TE_Orient = $strand;
            $desc = defined $desc ? "$desc " : '';
            my $TE_names = join( "/", @TEs );
            my $GFF_strand = $TE_Orient eq 'fwd' ? '+' : '-';
            if ($TE_Orient eq 'ref'){
              $GFF_strand = '.';
            }
            if ( $padding == 0 ) {
              my $header = "$TE_names $desc$ref:$start..$end $coords";
              if ( !$TE_Orient or $TE_Orient eq 'fwd' or $TE_Orient eq 'unk' or $TE_Orient eq 'ref' ) {
                print ">$name $TE_name.$TE_Orient $header\n$subseq\n";
                print GFF
                  "$name\t.\tchromosome\t$coords{flank_1_start}\t$coords{flank_2_end}\t.\t.\t.\tID=$name;Name=$name;Note=$header\n";
                print GFF
                  "$name\t.\ttransposable_element_insertion_site\t$coords{TE_start}\t$coords{TE_end}\t.\t$GFF_strand\t.\tID=$name.$start.$end;Name=$name;Note=$header\n";
              }
              if ( defined $revcomp_TE ) {
                if ( !$TE_Orient or $TE_Orient eq 'rc' or $TE_Orient eq 'unk' )
                {
                  print ">$name $TE_name.$TE_Orient $header\n$revcomp_subseq\n";
                  print GFF
                    "$name\t.\tchromosome\t$coords{flank_1_start}\t$coords{flank_2_end}\t.\t.\t.\tID=$name;Name=$name;Note=$header\n";
                  print GFF
                    "$name\t.\ttransposable_element_insertion_site\t$coords{TE_start}\t$coords{TE_end}\t.\t$GFF_strand\t.\tID=$name.$start.$end;Name=$name;Note=$header\n";
                }
              }
            } else {
              my $header = "$TE_names $desc$ref:$start..$end $coords";
              if ( !defined $name ) {
                $name = $header;
              }

              if ( !$TE_Orient or $TE_Orient eq 'fwd' or $TE_Orient eq 'unk' or $TE_Orient eq 'ref') {
                print ">$name $TE_name.$TE_Orient $header\n$subseq\n";
                print GFF
                  "$name\t.\tchromosome\t$coords{flank_1_start}\t$coords{flank_2_end}\t.\t.\t.\tID=$name;Name=$name;Note=$header\n";
                print GFF
                  "$name\t.\ttransposable_element_insertion_site\t$coords{TE_start}\t$coords{TE_end}\t.\t$GFF_strand\t.\tID=$name.$start.$end;Name=$name;Note=$header\n";
              }
              if ( defined $revcomp_TE ) {
                if ( !$TE_Orient or $TE_Orient eq 'rc' or $TE_Orient eq 'unk' )
                {
                  print ">$name $TE_name.$TE_Orient $header\n$revcomp_subseq\n";
                  print GFF
                    "$name\t.\tchromosome\t$coords{flank_1_start}\t$coords{flank_2_end}\t.\t.\t.\tID=$name;Name=$name;Note=$header\n";
                  print GFF
                    "$name\t.\ttransposable_element_insertion_site\t$coords{TE_start}\t$coords{TE_end}\t.\t$GFF_strand\t.\tID=$name.$start.$end;Name=$name;Note=$header\n";
                }
              }

            }

          }
      }
      else {
        warn "error retrieving $ref:$start..$end\n";
      }
    }
  }
}
}
### SOUBROUTINE ####
sub print_special_ping {
my ($start,$end,$seq,$name) = @_;
my $ref = 'Chr1';
my $before_ping = 'TACTCCCTCCAGTTCCTAATTTATTGGCGTTTTAGACATGATTACACATACCAACGAGTGATTAATTAGTGTGCCTACTTCCTGTTTTACCCCTAATAAATACGGTTTTGAGGTATACCAAATCATAACCTGCTCTTGGTTTGAGTGGCCATCCCATGCAAAAGTGAAGCTACCAGTTCTGAGTTTGATTTTTGGTGGCCACATTTTTTTGCATGGCACTGTGGACGCGTGACTGTTGGCATGCATGCACTTTTTTCAGGCGTTCACGCGCGCAGTCACCGGGAAAGATTTTTCACTCGCGCGCTGCTAA';
my $len_beforePing = length $before_ping;
my $after_ping = 'TAACAAAAATTTTGCCAGGTGCGTTTGCCGCCAATAACTCAAATTTTGGCGCTGGTTTGCAGCAGAATCCCCTTTGCCATCTCCATCTATTTATGCTGTGCAATACTAGCTAGTTTGACGCTCCAACTCCCTACATACTCCATCTGTTTAACATGAACTTAAACATGGAGCCATTACCTGCTCCATGGCGAAAATGGCGAAAGTAGTACAAGGAGACGAGTAGAAGGCAAAGAGCCAACAAGCAGGCATAGAAGAAGAAAATTACTCAGCGCGTGGGAGAGCAATGGTTGAACAGGAAACGAGAAGGCAAAGCATGCATGACGTTGGAAATAATGGGAGTAGTAATGGCCGGCTGCATGCATGCAGGGGCAAGCAGGTCCAAGATGCACCAAAAAAGCTAAAACGTCAAAAAATATCACACAATGCTTCAAGAGCTAAAACGTCAACTAATTTCAAACTGGAGGGAGTA';
my $stowaway_with_ping = 'TACTCCCTCCAGTTCCTAATTTATTGGCGTTTTAGACATGATTACACATACCAACGAGTGATTAATTAGTGTGCCTACTTCCTGTTTTACCCCTAATAAATACGGTTTTGAGGTATACCAAATCATAACCTGCTCTTGGTTTGAGTGGCCATCCCATGCAAAAGTGAAGCTACCAGTTCTGAGTTTGATTTTTGGTGGCCACATTTTTTTGCATGGCACTGTGGACGCGTGACTGTTGGCATGCATGCACTTTTTTCAGGCGTTCACGCGCGCAGTCACCGGGAAAGATTTTTCACTCGCGCGCTGCTAAGGCCAGTCACAATGGAGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATGATTATTTGCATGAAATGGGGATGAGAGAGAAGGAAAGAGTTTCATCCTGGTGAAACTCGTCAGCGTCGTTTCCAAGTCCTCGGTAACAGAGTGAAACCCCCGTTGAGGCCGATTCGTTTCATTCACCGGATCTCTTGCGTCCGCCTCCGCCGTGCGACCTCCGCATTCTCCCGCGCCGCGCCGCGCCACGCCTCCTTCCCGCGTGAACATTCCTCCTTCCCGCGCGAGCGATTCCACCATCTCCCCCGTCCGGCGCCTACGGAGTACACCGCAACCGGTCGCCCCAATCCGGCGCCTAGACCGTGACCCACCCGCCATCTTCCGCAAGACCGAATCCCCAACCCACCCACCATCTTCCGCCGCCCCCGTCCCCGTCCCCGGCCATGGATCCGTCGCCGGCCGTGGATCCGTCGCCGGCCGTGGATCCGTCGCCGGCTGCTGAAACCCGGCGGCGTGCAACCGGGAAAGGAGGCAAACAGCGCGGGGGCAAGCAACTAGGATTGAAGAGGCCGCCGCCGATTTCTGTCCCGGCCACCCCGCCTCCTGCTGCGACGTCTTCATCCCCTGCTGCGCCGACGGCCATCCCACCACGACCACCGCAATCTTCGCCGATTTTCGTCCCCGATTCGCCGAATCCGTCACCGGCTGCGCCGACCTCCTCTCTTGCTTCGGGGACATCGACGGCAAGGCCACCGCAACCACAAGGAGGAGGATGGGGACCAACATCGACCATTTCCCCAAACTTTGCATCTTTCTTTGGAAACCAACAAGACCCAAATTCATGGTACATGTATTTTCTTCTTTTTCTGTTACTTTCAACCTACGGTAACTCTAATTCATGGATGAGACTACTGCCATTGTGCAGTTCAATGCTTTTTCTTCATGTTATATTTCGTCCAGCTGTGAGTTATGGTTTGAAGATTGCTGTGGTTGTTTCATTGCTGAGTATGTGAAAGATAGATGGATGAAAGAGAGAATTATATTTTAGTCTGTAATCTTGCTCATCCAGTTGCTCATGTATGACCTTGGTTCTAGAATGTTGCCCTGACTGTATGCTTAATGTTCAGAGAAGTGATGCCTAAAGCAGTGAGATCAGTGGGATCAGATTAGCTATCGACATATAATATTAGCTATCTCAGTTGTGAAAGAGAGATGGGTGAAAAGGCACCCCTTGGATTAATTCTGTAGTATCAAATTCTGCACCTTGTCTGTCCATATGTTCTGCTTGGTTGGTGGGTGCAGTGCATTTGTAAAAAATAGTTTGCTTCTGATCCTTAATATATGTAACAGGGAATGAATTTTCACCCATCTCAGTTGTAAAGGTACTGTCTTGCTATGCAATATGTGTAAATTGACAAACCTGAAAATAGTCTGTTTGGAATTTGCAAAAGCAATTCGATAGTTTGGAATTTCCAAACCTCAGTCAGCAGTAGGCAATCCATTTTAGTTCTTGCTATGCACAAAAACAGTACACCTGATATGCTCATTTTAATACAACTTTTTTGTCTCTGTTACAGTTTGGTCAGGGGTTATCCTCCAGGAGGGTTTGTCAATTTTATTCAACAAAATTGTCCGCCGCAGCCACAACAGCAAGGTGAAAATTTTCATTTCGTTGGTCACAATATGGGATTCAACCCAATATCTCCACAGCCACCAAGTGCCTACGGAACACCAACACCCCAAGCTACGAACCAAGGCACTTCAACAAACATTATGATTGATGAAGAGGACAACAATGATGACAGTAGGGCAGCAAAGAAAAGATGGACTCATGAAGAGGAAGAGAGACTGGTATTCATCGGATACTTTTACATTTCCATATGTCTTTGTTTTGACTAATACTTGACAGGTCATTAACTGATTCTTGTAGGCCAGTGCTTGGTTGAATGCTTCTAAAGACTCAATTCATGGGAATGATAAGAAAGGTGATACATTTTGGAAGGAAGTCACTGATGAATTTAACAAGAAAGGGAATGGAAAACGTAGGAGGGAAATTAACCAACTGAAGGTTCACTGGTCAAGGTTGAAGTCAGCGATCTCTGAGTTCAATGACTATTGGAGTACGGTTACTCAAATGCATACAAGCGGATACTCCGACGACATGCTTGAGAAAGAGGCACAGAGGCTGTATGCAAACAGGTTTGGAAAACCTTTTGCGTTGGTCCATTGGTGGAAGATACTCAAAGATGAGCCCAAATGGTGTGCTCAGTTTGAATCAGAGAAAGACAAGAGCGAAATGGATGCTGTTCCAGAACAGCAGTCACGTCCTATTGGTAGAGAAGCAGCAAAGTCTGAGCGCAATGGAAAGCGCAAGAAAGAAAATGTTATGGAAGGCATTGTCCTCCTAGGGGACAATGTCCAGAAAATTATAAAGGTCCACGAAGACCGGAGGGTGGATCGTGAAAAGGCCACCGAAGCACAGATTCAGATATCAAATGCAACATTGTTGGCCGCTAAGGAGCAGAAGGAAGCAAAGATGTTCGATGTGTACAATACTCTATTAAGTAAGGATACAAGCAACATGTCTGAAGATCAAATGGCTAGCCACCAGAGGGCAATACGGAAATTAGAGGAGAAGCTATTTGCGGATTAAGGTGAGTTTTATAAACTGACCACTATTTTCTGAAATGTATGAATTCTGAAATTTATATACAATTGTGTAAACATGGAAAATTAGATAATGTATGCATGATGCACAACATGTGCGTGCAGCACTATTTAATGGCAGTTTCACAAGTGTGAAAACTGACCACTATAGTACTATTGTGGTGTGAAAACTGACCACTACTATTGTGGTGTGAATGCTACTGTGGTGTGAAAACTGACCACTATAGTTTCACATTCCTGGATGCAGCCCTCCTCTATATATATAGATACAGTCCTCATCTCTTCCTGGCATACACACAGCCCTCTTCTCTAATTCCTGGACGCAGTCCTCATCTCTTCCTGGCATAGACGCAGCCCTTCTCTCTTCCTGTTTAGTTCAACAACATTGAGGTGATCTGCCTTTCTTTGAAGTTTCTATCTTTTTTCACTGCTGTGAATGATTATTTCTCTGCTGTGAATGATTATTTCTCCAATCTTCCTTTGTTCACCTTCTCTCTTTCTCTGCTGTGAAGATGTCTGGAAATGAAAATCAGATTCCTGTGTCCTTGTTGGACGAGTTTCTCGCTGAGGATGAGATCATGGATGAGATAATGGATGATGTTCTCCATGAAATGATGGTGTTATTGCAGTCCTCCATCGGAGATCTTGAAAGAGAGGCTGCTGACCATCGTTTGCATCCAAGGAAGCACATCAAGAGGCCACGAGAGGAAGCACATCAAAATTTGGTGAATGATTATTTCTCTGAAAATCCTCTATATCCTTCCAATATTTTTCGCCGAAGATTTCGTATGTACAGGCCGCTGTTTTTACGTATTGTGGACGCATTAGGCCAGTGGTCAGATTACTTTACTCAGAGGGTAGATGCCGCTGGTAGGCAAGGGCTTAGTCCATTACAAAAGTGTACTGCAGCAATTCGCCAATTGGCTACTGGTAGTGGTGCTGATGAACTAGATGAGTATTTGAAGATTGGAGAGACTACTGCTATGGATGCTATGAAAAATTTTGTGAAAGGAATTAGAGAAGTATTTGGTGAAAGATATCTCAGGCGTCCCACTGTAGAAGATACTGAACGACTACTCGAGCTTGGTGAGAGACGCGGTTTTCCTGGTATGTTCGGTAGCATTGACTGTATGCATTGGCAATGGGAAAGGTGCCCAACTGCGTGGAAGGGTCAGTTCACTCGTGGTGATCAAAAAGTGCCAACGCTGATTCTTGAGGCAGTGGCATCACATGATCTTTGGATTTGGCATGCGTTCTTTGGAGTAGCAGGTTCTAACAATGATATCAATGTTTTGAGCCGATCTACTGTGTTTATCAATGAGCTGAAAGGACAAGCTCCTAGAGTGCAGTACATGGTAAATGGGAATCAATACAACGAAGGTTATTTTCTTGCTGATGGAATTTACCCTGAATGGAAGGTATTTGCTAAGTCATATCGACTCCCTATCACTGAGAAGGAGAAGTTGTATGCACAACATCAAGAAGGGGCAAGAAAGGATATCGAGAGAGCATTTGGTGTTCTACAACGTCGATTCTGCATCTTAAAACGACCAGCCCGTCTATATGACCGAGGTGTACTCCGTGATGTTGTCCTAGGTTGCATCATACTTCACAATATGATAGTTGAAGATGAGAAGGAAGCGCGACTTATTGAAGAAAATCTAGATTTAAATGAGCCTGCTAGTTCATCAACGGTTCAGGCACCAGAATTCTCTCCTGACCAGCATGTTCCATTAGAAAGAATTTTAGAAAAGGATACTAGTATGAGAGATCGTTTGGCTCATCGCCGACTCAAGAATGATTTGGTGGAACATATATGGAATAAGTTTGGTGGTGGTGCACATTCATCTGGTAATTATGTTTTTATTTTGCATTATTAGTTATCTATGGTACTAAGATATGTACAAGTTTCTCTAAATTGCACTAAATCTGTGGTTCATATTGGATATGTGTAAACTATGAATGTAGCCTGACTAAAACCATCATTCATGCTGAACTGGTTTTTGTTTTGTATATGCAGGATGAAACAAGGAACTAGGTTTCTGAACGCATTACGGACTGAAGGTTGAGGGGCAGAATGATCCACCCAGTTGCTTCTATCAGATCACTAAAGTTTCATTTCACTGTTTTATTTTGGACACTTGATGCTTGTGTGCATCCGATGAATGTTTAATTTGGTCACCTGATGCTTGTGTGCATCCGATGAATGTTTAATTTGGTCACCTGATGCTTGTATGCAGTTATCTATCTTATTTGTTAATGTTGCTGGTACTGAGGATTTTTAGAAGTGAAATGCACAAGTTGCTGTGTTTTTTGACTGATCCTTGTGTGCACTTGACGTTGTATGTGACAAATGATGGTTCCCAGTTGTGCACCTGATTCATGATTCAGTTATTCAGTTTAAATTGACGTTGTTTGTGTGCACCTTTTGTCAGTTAGCCAGTTACGGCTGGAAGTTGTGTAAGTTTGTGTGACGCCTGGCTACAGGATTTTGGGTACAAATGATCCCAGCAACTTGTATCAATTAAATGCTTTGCTTAGTCTTGGAAACGTCAAAGTGAAACCCCTCCACTGTGGGGATTGTTTCATAAAAGATTTCATTTGAGAGAAGATGGTATAATATTTTGGGTAGCCGTGCAATGACACTAGCCATTGTGACTGGCCTAACAAAAATTTTGCCAGGTGCGTTTGCCGCCAATAACTCAAATTTTGGCGCTGGTTTGCAGCAGAATCCCCTTTGCCATCTCCATCTATTTATGCTGTGCAATACTAGCTAGTTTGACGCTCCAACTCCCTACATACTCCATCTGTTTAACATGAACTTAAACATGGAGCCATTACCTGCTCCATGGCGAAAATGGCGAAAGTAGTACAAGGAGACGAGTAGAAGGCAAAGAGCCAACAAGCAGGCATAGAAGAAGAAAATTACTCAGCGCGTGGGAGAGCAATGGTTGAACAGGAAACGAGAAGGCAAAGCATGCATGACGTTGGAAATAATGGGAGTAGTAATGGCCGGCTGCATGCATGCAGGGGCAAGCAGGTCCAAGATGCACCAAAAAAGCTAAAACGTCAAAAAATATCACACAATGCTTCAAGAGCTAAAACGTCAACTAATTTCAAACTGGAGGGAGTA';my $len_afterPing = length $after_ping;
#my $pre_TE_len = $padding - $len_beforePing + 3;
#my $post_TE_len = $padding - $len_afterPing + 3;
my $pre_TE_len = $padding - $len_beforePing ;
my $post_TE_len = $padding - $len_afterPing ;

my $pre_TE = substr $seq, ($start - ($pre_TE_len)-1) ,$pre_TE_len;
my $post_TE = substr $seq, $end ,$post_TE_len;
my $special_ping_plus_flank = $pre_TE . $stowaway_with_ping . $post_TE;
my $total_len = length $special_ping_plus_flank;
my $stowaway_start = $pre_TE_len + 1;
my $stowaway_break_end = $padding;
my $ping_start = $padding + 1;
my $ping_end = $total_len - $padding;#length $pre_TE . $stowaway_with_ping;
my $stowaway_break_start = $ping_end + 1;
my $stowaway_end = $ping_end + $len_afterPing;
my $flank2_start = $ping_end +1;
my $info = "$ref:$start..$end FLANK1:1..$padding STOWAWAY_prePing:$stowaway_start..$stowaway_break_end TE:$ping_start..$ping_end STOWAWAY_postPing:$stowaway_break_start..$stowaway_end FLANK2:$flank2_start..$total_len";
print ">$name HEG4/EG4/A119/A123 $info\n$special_ping_plus_flank\n";
print GFF "$name\t.\tchromosome\t1\t$total_len\t.\t.\t.\tID=$name;Name=$name;Note=shared:HEG4/EG4/A119/A123 $info\n";
print GFF "$name\t.\ttransposable_element_insertion_site\t$ping_start\t$ping_end\t.\t.\t.\tID=$name.$start.$end;Name=$name;Note=ping HEG4/EG4/A119/A123 $info\n";
}
# my $header = "$TE_names $desc$ref:$start..$end $coords";
#"$name\t.\ttransposable_element_insertion_site\t$coords{TE_start}\t$coords{TE_end}\t.\t$GFF_strand\t.\tID=$name.$start.$end;Name=$name;Note=$header\n";
