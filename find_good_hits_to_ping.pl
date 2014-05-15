#!/usr/bin/perl -w
use strict;

my $sam = shift;
my $flank = 600;

my @sam_path = split /\// , $sam;
my $SAM = pop @sam_path;
my ($lib) = $SAM =~ /(.+)\.sam/;

open INSAM, $sam or die "Can't open $sam\n";
open GOOD , ">goodTEHits.$lib.txt" or die "Can't open goodTEHits.$lib.txt for writing\n";
open READS , ">goodTEHits.$lib.sam" or die "Can't open goodTEHits.$lib.sam for writing\n";
#open BAD , ">notTEHits.txt" or die "Can't open notTEHits.txt for writing\n";

my %pairs;
my %refs;
print READS "\@CO sample:$lib\n";
while (my $line = <INSAM>){
  chomp $line;
  if ($line =~ /\@SQ\s|\@HD\s|\@RG\s|\@PG\s/){
    next unless $line =~ /^\@SQ\sSN/;
    # @SQ     SN:ping5:HEG4/EG4       LN:6547
    my ($ref,$len) = $line =~ /\@SQ\s+SN:(\S+)\s+LN:(\d+)/;
    $refs{$ref}=$len;
    print READS "$line\n";
    next;
  }
  if ($line =~ /XM:i:(\d+)/){
    next if $1 > 3;
  }
  my @line = split /\t/ , $line;
  next if $line[2] eq '*';
  next if $line[2] !~ /ping/;
  my $read_name = $line[0];
  my $read_pair;
  if ($line[1] & 64){
    $read_pair = 'first';
  }else{
    $read_pair = 'second';
  }
  $pairs{$line[0]}{$read_pair}{RNAME}=$line[2];
  $pairs{$line[0]}{$read_pair}{POS}=$line[3];
  $pairs{$line[0]}{$read_pair}{MAPQ}=$line[4];
  $pairs{$line[0]}{$read_pair}{CIGAR}=$line[5];
  $pairs{$line[0]}{$read_pair}{read_len}=length $line[9];
  $pairs{$line[0]}{$read_pair}{LINE}=$line;
}
#my $ping_len = 5341; 
my %ping_count;
my @bothInside;
my $uniq_ping_s = 252; ## mping and ping same upto 251 of ping
my $uniq_ping_e = 180; ## mping and ping last 180 of ping
foreach my $read (keys %pairs){
  next if scalar keys %{$pairs{$read}} != 2;
  ## check to see if both pairs map to same ping
  my $ref_first = $pairs{$read}{first}{RNAME};
  my $ref_second = $pairs{$read}{second}{RNAME};
  next if $ref_first ne $ref_second;
  
  ## one needs to be outside the ping,
  ## one needs to be inside the ping
  ## ---flank------(ping)----flank---
  my $start_first = $pairs{$read}{first}{POS}; 
  my $start_second = $pairs{$read}{second}{POS}; 
  my $len_first = $pairs{$read}{first}{read_len}; 
  my $len_second = $pairs{$read}{second}{read_len}; 
  my $cigar_first = $pairs{$read}{first}{CIGAR}; 
  my $cigar_second = $pairs{$read}{second}{CIGAR}; 
  my @lens = ($len_first,$len_second);
  my @cigars = ($cigar_first,$cigar_second);
  if ($start_first > $start_second){
    @lens= ($len_second,$len_first);
    @cigars= ($cigar_second,$cigar_first);
  }
  my @starts = sort {$a <=> $b}($start_first,$start_second);
  my $ref_len = $refs{$ref_first};
  my $ref_start = $flank+1;
  my $ref_end = $ref_len - $flank;
  my $uniq_ref_start = $flank + $uniq_ping_s + 1;
  my $uniq_ref_end = $ref_len - $flank - $uniq_ping_e;
  
  #my $insertSize = $starts[1] - ($starts[0]+$lens[0]);
  my $insertSize = ($starts[1]+$lens[1]) - $starts[0];
  my $in_uniq = 'notUniq';
  if ($starts[0] <= $ref_start - $lens[0] and ( $starts[1] >= $ref_start and $starts[1] <= ($ref_end-$lens[1]))){
  ## if first is before ping, second should be inside
    #print "first before, second in\n";
    $in_uniq = 'inUniq' if ($starts[1] >= $uniq_ref_start and $starts[1] <= $uniq_ref_end-$lens[1]);
    print GOOD "$sam\t$ref_first\tfirstBefore,secondIn($in_uniq)\tis:$insertSize\n";
    print READS "\@CO firstBefore, secondIn($in_uniq): $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    print READS $pairs{$read}{first}{LINE},"\n";
    print READS $pairs{$read}{second}{LINE},"\n";
    $ping_count{$ref_first}++;
  }elsif(( $starts[0] >= $ref_start and $starts[0] <= ($ref_end - $lens[0])) and $starts[1] >= $ref_end ) {
  ## if first is inside ping, second should be after
    #print "first in, second after\n";
    $in_uniq = 'inUniq' if ($starts[0] >= $uniq_ref_start and $starts[0] <= $uniq_ref_end-$lens[0]);
    print GOOD "$sam\t$ref_first\tfirstIn($in_uniq),secondAfter\tis:$insertSize\n";
    print READS "\@CO firstIN, secondAfter($in_uniq): $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    print READS $pairs{$read}{first}{LINE},"\n";
    print READS $pairs{$read}{second}{LINE},"\n";
    $ping_count{$ref_first}++;

  }elsif($starts[0] >= $uniq_ref_start and $starts[0] <= ($uniq_ref_end-$lens[0]) and $starts[1] >= $uniq_ref_start and $starts[1] <= $uniq_ref_end-$lens[1]){
    next if $insertSize < 100; 
    warn "BothInsidePing: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end,uniqToPing:$uniq_ref_start..$uniq_ref_end\n";
    #if ($lib =~ /NB/ and $ref_first =~ /NB/){
      warn "$sam\t$ref_first\tbothInsidePing\tis:$insertSize\n";
      push @bothInside, "\@CO BothInsidePing: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end,uniqToPing:$uniq_ref_start..$uniq_ref_end\n";
      push @bothInside, $pairs{$read}{first}{LINE},"\n";
      push @bothInside,  $pairs{$read}{second}{LINE},"\n";
    #}
  }
  elsif ($starts[0] >= $ref_end and $starts[1] >= $ref_end){
    warn "BothAfter: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    #print READS "\@CO BothAfter: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    #print READS $pairs{$read}{first}{LINE},"\n";
    #print READS $pairs{$read}{second}{LINE},"\n\n";
  }
  elsif ($starts[0] <= $ref_start-$lens[0] and $starts[1]<=$ref_start-$lens[1]){
    warn "BothBefore: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    #print READS "\@CO BothBefore: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    #print READS $pairs{$read}{first}{LINE},"\n";
    #print READS $pairs{$read}{second}{LINE},"\n\n";
  }
  elsif (  (($starts[0] < ($ref_start-30)) and ($starts[0]+$lens[0]) > ($ref_start+30)   ) and ($starts[1]>$ref_start and $starts[1] <= $ref_end) ){
    ## firstJunction, secondIn
    my @cigar_M0 = $cigars[0]  =~ /(\d+)M/g; 
    my $cigar_M0;
    foreach my $M (@cigar_M0){
      $cigar_M0+=$M;
    }
    my @cigar_M1 = $cigars[1]  =~ /(\d+)M/g; 
    my $cigar_M1;
    foreach my $M (@cigar_M1){
      $cigar_M1+=$M;
    }
    if ($lens[0]*.9 < $cigar_M0 and $lens[1]*.9 < $cigar_M1){
      print GOOD "$sam\t$ref_first\tfirstJunction,secondIn\n";
      print READS "firstJuctiion, secondIn: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
      print $pairs{$read}{first}{LINE},"\n";
      print $pairs{$read}{second}{LINE},"\n\n";

    }else{
      warn "Neither:almostJunctionCase1: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    }
  }
  elsif (  ($starts[0] >= $ref_start and $starts[0] <= $ref_end-$lens[0]) and ($starts[1] <= ($ref_end-30) and ($starts[1]+$lens[1]) >= ($ref_end+30))){
    ## firstIn, secondJunction
    my @cigar_M0 = $cigars[0]  =~ /(\d+)M/g; 
    my $cigar_M0;
    foreach my $M (@cigar_M0){
      $cigar_M0+=$M;
    }
    my @cigar_M1 = $cigars[1]  =~ /(\d+)M/g; 
    my $cigar_M1;
    foreach my $M (@cigar_M1){
      $cigar_M1+=$M;
    }
    if ($lens[0]*.9 < $cigar_M0 and $lens[1]*.9 < $cigar_M1){
      print GOOD "$sam\t$ref_first\tfirstIn,secondJunction\n";
      print READS "firstIN, secondJunction: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
      #print warn "firstIN, secondJunction: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
      print $pairs{$read}{first}{LINE},"\n";
      print $pairs{$read}{second}{LINE},"\n\n";
    }else{
      warn "Neither:almostJunctionCase2: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    }
  }
  # firstSecondJunction, secondAfter
  # -------[      ]-----
  #             ---- -----
  elsif( ($starts[0] <= ($ref_end - 30) and ($starts[0]+$lens[0] > $ref_end+30)) and ($starts[1] > $ref_end )){
    my @cigar_M0 = $cigars[0]  =~ /(\d+)M/g;
    my $cigar_M0;
    foreach my $M (@cigar_M0){
      $cigar_M0+=$M;
    }
    my @cigar_M1 = $cigars[1]  =~ /(\d+)M/g;
    my $cigar_M1;
    foreach my $M (@cigar_M1){
      $cigar_M1+=$M;
    }
    if ($lens[0]*.9 < $cigar_M0 and $lens[1]*.9 < $cigar_M1){
      print GOOD "$sam\t$ref_first\tfirstIn,secondJunction\n";
      #print warn "firstSecondJunction, secondAfter: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
      print READS "firstSecondJunction, secondAfter: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
      print $pairs{$read}{first}{LINE},"\n";
      print $pairs{$read}{second}{LINE},"\n\n";
    }else{
      warn "Neither:almostJunctionCase3: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    }
  }
  # firstBefore, secondJunction
  # -------[      ]-----
  # ---- -----
  elsif ($starts[0] <= ($ref_start - $lens[0]) and ($starts[1] >= ($ref_start-30)) and (($starts[1]+$lens[1]) >= ($ref_start+30))  ){
    my @cigar_M0 = $cigars[0]  =~ /(\d+)M/g;
    my $cigar_M0;
    foreach my $M (@cigar_M0){
      $cigar_M0+=$M;
    }
    my @cigar_M1 = $cigars[1]  =~ /(\d+)M/g;
    my $cigar_M1;
    foreach my $M (@cigar_M1){
      $cigar_M1+=$M;
    }
    if ($lens[0]*.9 < $cigar_M0 and $lens[1]*.9 < $cigar_M1){
      print GOOD "$sam\t$ref_first\tfirstIn,secondJunction\n";
      #print warn "firstBefore, secondJunction: $starts[0]($lens[0]),$starts[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
      print READS "firstBefore, secondJunction: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),is:$insertSize,$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
      print $pairs{$read}{first}{LINE},"\n";
      print $pairs{$read}{second}{LINE},"\n\n";
    }else{
      warn "Neither:almostJunctionCase4: $starts[0]:$cigars[0]($lens[0]),$starts[1]:$cigars[1]($lens[1]),$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
    }
  }
 
  else {
    warn "Neither: $starts[0]($lens[0]),$starts[1]($lens[1]),$ref_first($refs{$ref_first}):$ref_start..$ref_end\n";
  }
}
if ((keys %ping_count) == 1){
  print join ("\n", @bothInside);
}
