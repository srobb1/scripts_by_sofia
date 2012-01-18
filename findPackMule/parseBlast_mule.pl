#!/usr/bin/perl -w
use strict;
use Bio::SearchIO; 
use Bio::Range;
my $blast = shift; #blastout of db_file against packmule TSD-TIR-TIR-TSD
my $db_file = shift; #same file used as a database in the above blast
my $db2 = shift; # nr, species specific of the above organsims
my $min_insert = shift ;#2000;
my $max_insert = shift ;#8000;
my $min_ident  = shift ;#85;
my $query_half = shift ;#0.4;

my $searchIO_obj = new Bio::SearchIO(-format => 'blast', 
                           -file   => $blast);

open OUTMATCH , ">match.table.txt";
my %inserts;
print "h_name\tper_iden\thsp_len\th_strand\th_start\th_end\tq_start\tq_end\n";
while( my $result = $searchIO_obj->next_result ) {
  my $q_len = $result->query_length;
  while( my $hit = $result->next_hit ) {
    my %hash;
    my $hsp_count = 0;
    my $biggest_end = 0;
    while( my $hsp = $hit->next_hsp ) {
      if( $hsp->length('total') > $q_len * $query_half ) {
        if ( $hsp->percent_identity >= $min_ident ) {
		my $h_end =  $hsp->end('hit');
		$biggest_end =  $h_end > $biggest_end ? $h_end : $biggest_end;
        }
      }
    }
    $hit->rewind;
    while( my $hsp = $hit->next_hsp ) {
      my $add = 1;
      if( $hsp->length('total') > $q_len * $query_half ) {
        if ( $hsp->percent_identity >= $min_ident ) {
          my $q_start = $hsp->start('query');
          my $q_end = $hsp->end('query');
          my $h_start = $hsp->start('hit');
          my $h_end = $hsp->end('hit');
          my $h_name = $hit->name;
          my $h_len = $hit->length;
          my $per_iden = $hsp->percent_identity;
          my $hsp_len = $hsp->length('total');
          my $h_strand = $hsp->strand('hit');
          print "$h_name\t$per_iden\t$hsp_len\t$h_strand\t$h_start\t$h_end\t$q_start\t$q_end\n";
          if (!exists $hash{$h_name}{$h_strand}){
            for (my $i = 0 ; $i < $biggest_end+1 ; $i++){
              $hash{$h_name}{$h_strand}[$i] = '-';
            }
          }
          my $j;
          if ( abs($q_end - $q_len*.5) <= ($q_len*.5 - ($q_len * $query_half))){
            $j = 1;
          }else{
            $j = 2;
          }
          for (my $i = $h_start ; $i < $h_end+1 ; $i++){
            $hash{$h_name}{$h_strand}[$i] = $j;
          }
        }
      }
    }
    foreach my $h_name (sort keys %hash){  
    foreach my $strand (keys %{$hash{$h_name}}){
      my $hits = join '' , @{$hash{$h_name}{$strand}};
      my ($first,$second);
      if ($strand > 0){
        $first = 1;
        $second = 2;
      } else {
        $first = 2;
        $second = 1;
      }
      my $pattern =  '(('.$first.'{' . $q_len*$query_half . ',})(\D{' . $min_insert . ',' . $max_insert . '})('.$second.'{' . $q_len*$query_half . ',}))';
      #$hits =~ s/\D+$//;
      while ($hits =~ /$pattern/g){
        my $len_TIR1 = length $2;
        my $len_TIR2 = length $4;
        my $insert_len = length $3;
        my $total_len = length $1;
        my $end = pos ($hits) - 1;
        my $start = $end - $total_len +1 ; 
        my $insert_start =  $start + $len_TIR1;
        my $insert_end =    $insert_start + $insert_len - 1;
        my $range = Bio::Range->new(-start=>$start, -end=>$end, -strand=>$strand);
        my $add = 1;
        foreach my $id (keys %{$inserts{$h_name}}){
            my $stored_range  = $inserts{$h_name}{$id}{range};
            if ($range->overlaps($stored_range)){
              if ($range->length > $stored_range->length){
                ## new range is longer, throw out old one
                delete  $inserts{$h_name}{$id};
                #add it: $ranges{$h_name}{"$start..$end"}=$range;
              }elsif ($range->length < $stored_range->length){
                ## old range is bigger,
                ## don't add the new one
                $add = 0;
              }elsif (($range->length == $stored_range->length)){
                ##don't add it
                $add = 0;
              }
            }
            else {
            ##add it: $ranges{$h_name}{"$start..$end"}=$range;
            }
         
        }
        if ($add){
            $inserts{$h_name}{"$start..$end"}{range}=$range;
            $inserts{$h_name}{"$start..$end"}{line}="$h_name\t$strand\t$total_len\t$start\t$insert_start\t$insert_end\t$end\n";
            $inserts{$h_name}{"$start..$end"}{t_start}=$start;
            $inserts{$h_name}{"$start..$end"}{i_start}=$insert_start;
            $inserts{$h_name}{"$start..$end"}{i_end}=$insert_end;
            $inserts{$h_name}{"$start..$end"}{t_end}=$end;
            $inserts{$h_name}{"$start..$end"}{t_len}=$total_len;
            $inserts{$h_name}{"$start..$end"}{i_len}=$insert_len;
        }
        #print OUTMATCH "$h_name\t$strand\t$total_len\t$start\t$insert_start\t$insert_end\t$end\n";
      }
    }
    }
  }
}

print "\n\n";
my $fasta;
my $insert;
my $file = $db_file;
open OUTPACK , ">packmule.fa";
open OUTINSERT , ">insertOnly.fa";
print OUTMATCH "ref\tstrand\ttotal_len\tstart\tinsert_start\tinsert_end\tend\n";
print "\n\nPackMule inserts:\n";
print "ref\tstrand\ttotal_len\tstart\tinsert_start\tinsert_end\tend\n";

foreach my $h_name (sort keys %inserts){
  my @ranges = keys %{$inserts{$h_name}} ;
  foreach my $range (keys %{$inserts{$h_name}}){
    my $start =$inserts{$h_name}{$range}{t_start};
    my $i_start=$inserts{$h_name}{$range}{i_start};
    my $i_end = $inserts{$h_name}{$range}{i_end};
    my $end = $inserts{$h_name}{$range}{t_end};
    my $t_len = $inserts{$h_name}{$range}{t_len};
    my $line = $inserts{$h_name}{$range}{line};
    print $line;
    print OUTMATCH $line;
    my $ref = `fastacmd -s "$h_name" -d $db_file`;
    my ($header,@seq) = split /\n/, $ref;
    my $ref_seq = join '',@seq;
    my $insert_len = $i_end - $i_start +1 ;
    my $packmule_seq = substr($ref_seq, $start-1,$t_len);
    my $insert_seq = substr($ref_seq, $i_start-1, $insert_len);
    print OUTPACK  ">$h_name:$start..$end insert:$i_start..$i_end total_len=$t_len insert_len=$insert_len\n$packmule_seq\n";
    print OUTINSERT ">$h_name:$i_start..$i_end len=$insert_len\n$insert_seq\n";
    $fasta .=  ">$h_name:$start..$end insert:$i_start..$i_end total_len=$t_len insert_len=$insert_len\n$packmule_seq\n";
    $insert .= ">$h_name:$i_start..$i_end len=$insert_len\n$insert_seq\n";
  }
}
print "\n\n##FASTA:Packmule\n$fasta";
print "\n\n##FASTA:insert\n$insert";

