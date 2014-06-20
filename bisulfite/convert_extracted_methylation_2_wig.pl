#!/usr/bin/perl -w
use strict;

#my $dir = shift;

##head CpG_OB_RIL43.RIL43_p1.fq_bismark_bt2_pe.txt
#Bismark methylation extractor version v0.9.0
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	376063	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	376060	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375985	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375911	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375914	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1879:2169_1:N:0:CTTGTA	-	Chr8	375985	z
#DGGXHXP1:377:C2K44ACXX:8:1101:1903:2170_1:N:0:CTTGTA	+	Chr1	37611389	Z
#DGGXHXP1:377:C2K44ACXX:8:1101:1903:2170_1:N:0:CTTGTA	+	Chr1	37611384	Z
#DGGXHXP1:377:C2K44ACXX:8:1101:1903:2170_1:N:0:CTTGTA	+	Chr1	37611343	Z

my $file = shift;    #<$dir/C*txt>;
my $fasta = shift;
my $i    = 0;
my %methy;
my @path = split /\//, $file;
my $file_name = pop @path;

open FASTA, $fasta or  die "Can't open $fasta";
my %seq;
my $id;
while ( my $line = <FASTA> ) {
  chomp $line;
  if ($line =~ /^>(\S+)/){
    $id = $1;
  }elsif($line !~ /^\s*$/) {
    $line =~ s/\s+//g;
    push @{$seq{ucfirst($id)}{seq}} ,  (split ('', $line));
  }
}
foreach my $id (keys %seq){
  my $seq = join '' , @{$seq{$id}{seq}};
  my @triplets_1 = $seq =~ /(.{3})/g;
  my @triplets_2 = (substr ($seq, 1)) =~ /(.{3})/g;
  my @triplets_3 = (substr ($seq, 2)) =~ /(.{3})/g;
  my @di_1 = $seq =~ /(.{2})/g;
  my @di_2 = (substr ($seq, 1))=~ /(.{2})/g;
  for (my $i = 0 ; $i < @triplets_1 ; $i++){
    my $trip = $triplets_1[$i];
    my $pos = ($i * 3) + 1;
    if($trip =~ /C[^G][^G]/){
      $seq{$id}{met}{$pos}='CHH';
    }elsif ($trip =~ /C[^G]G/){
      $seq{$id}{met}{$pos}='CHG';
    }
  }
  for (my $i = 0 ; $i < @triplets_2 ; $i++){
    my $trip = $triplets_2[$i];
    my $pos = ($i * 3) + 2;
    if($trip =~ /C[^G][^G]/){
      $seq{$id}{met}{$pos}='CHH';
    }elsif ($trip =~ /C[^G]G/){
      $seq{$id}{met}{$pos}='CHG';
    }
  }
  for (my $i = 0 ; $i < @triplets_3 ; $i++){
    my $trip = $triplets_3[$i];
    my $pos = ($i * 3) + 3;
    if($trip =~ /C[^G][^G]/){
      $seq{$id}{met}{$pos}='CHH';
    }elsif ($trip =~ /C[^G]G/){
      $seq{$id}{met}{$pos}='CHG';
    }
  }
  for (my $i = 0 ; $i < @di_1 ; $i++){
    my $di = $di_1[$i];
    my $pos = ($i * 2) + 1;
    if($di =~ /CG/){
      $seq{$id}{met}{$pos}='CpG';
    }
  }
  for (my $i = 0 ; $i < @di_2 ; $i++){
    my $di = $di_2[$i];
    my $pos = ($i * 2) + 2;
    if($di =~ /CG/){
      $seq{$id}{met}{$pos}='CpG';
    }
  }
}


#CHH_OB_RIL16_ping.RIL16_ping_p1.fq_bismark_bt2_pe.txt
#CHG_RIL16_ping.bismarkextract.out
my ( $type, $sample );
if (
    $file_name =~ /^(CpG|CHH|CHG)_(?:OT|OB)_(.+)\.\2_p1.fq_bismark_bt2_pe.txt/ )
{
  $type   = $1;
  $sample = $2;
} elsif ( $file_name =~ /^(CpG|CHH|CHG)_(.+)\.bismarkextract.out/ ) {
  $type   = $1;
  $sample = $2;
}
if ( !-e "$file.sorted" and $file !~ /sorted$/ ) {
  system("sort -k3 $file | grep -v Bismark > $file.sorted");
  $file = "$file.sorted";
}
print
  "track type=wiggle_0 name=\"$sample.$type\" description=\"$type in $sample\"\n";
open IN, $file or die "Can't open $file for reading $!\n";
my $last = '';
while ( my $line = <IN> ) {
  chomp $line;
  next if $line =~ /^Bismark methylation extractor version/;
  my ( $read, $strand, $ref, $pos, $methy ) = split /\t/, $line;
  $strand = 0;
  $ref = ucfirst $ref if $ref =~ /^Chr/i;
  $methy{$sample}{$ref}{$pos}{$type}{$methy}++;
  if ( ( $ref ne $last and $last ne '' ) or eof(IN) ) {

    #foreach my $sample ( keys %methy ) {
      #foreach my $ref ( keys %{ $methy{$sample} } ) {
        foreach my $pos ( sort { $a <=> $b } keys %{ $methy{$sample}{$last} } ) {
          $i++;

          foreach my $type ( keys %{ $methy{$sample}{$last}{$pos} } ) {
            my @methy_codes =
              keys %{ $methy{$sample}{$last}{$pos}{$type} };

            my $code   = $methy_codes[0];
            my $big    = uc $code;
            my $little = lc $code;
            my ( $methy_yes, $methy_no ) = ( 0, 0 );
            if ( exists $methy{$sample}{$last}{$pos}{$type}{$big} ) {
              $methy_yes = $methy{$sample}{$last}{$pos}{$type}{$big};
            }
            if ( exists $methy{$sample}{$last}{$pos}{$type}{$little} ) {
              $methy_no = $methy{$sample}{$last}{$pos}{$type}{$little};
            }

            my $total = $methy_yes + $methy_no;
            my $percent = ( $methy_yes / ($total) ) * 100;
            if ( $percent > 0 and $percent < 0.1 ) {
              $percent = 0.1;
            }
            if ( $percent == 0 ) {
              $percent = -100;
            }
            my $pretty_percent = sprintf( '%.1f', $percent );
            my $s = $pos - 1;
            my $nt = uc(${$seq{ucfirst $ref}{seq}}[$s]);
            my $met_type = $seq{ucfirst $ref}{met}{$pos};
            print "$last\t$s\t$pos\t$pretty_percent\n" if uc $nt eq 'C' and $met_type eq $type;;

          }
        }
        delete $methy{$sample}{$last};
        #print "Deleting $last\n"; 
        $last = $ref;
      #}
    #}
  } elsif ( $last eq '' ) {
    $last = $ref;
  }
}
