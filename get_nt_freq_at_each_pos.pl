#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;

my $seqIO_obj = Bio::SeqIO->new (-file => $file , -format => 'fasta');
my %nt_freq;
while (my $seq_obj = $seqIO_obj->next_seq){
	my $seq = $seq_obj->seq;
	$seq = uc $seq;
	my @seq = split '' , $seq;
	for (my $i = 0 ; $i < @seq ; $i++){
		$nt_freq{$i}{$seq[$i]}++;
	}
}

print "pos\tdepth\tfreq_A\tfreq_T\tfreq_G\tfreq_C\tfreq_N\n";
foreach my $pos (sort {$a <=> $b} keys %nt_freq){
	my $total_count = 0;
	foreach my $nt (sort keys %{$nt_freq{$pos}}){
		my $count = $nt_freq{$pos}{$nt};
		$total_count += $count;
	}
		my $count_A = $nt_freq{$pos}{A};
		my $count_T = $nt_freq{$pos}{T};
		my $count_G = $nt_freq{$pos}{G};
		my $count_C = $nt_freq{$pos}{C};
		my $count_N = $nt_freq{$pos}{N};
		my $freq_A = defined $count_A ? $count_A / $total_count : 0;
		my $freq_T = defined $count_T ? $count_T / $total_count : 0;
		my $freq_G = defined $count_G ? $count_G / $total_count : 0;
		my $freq_C = defined $count_C ? $count_C / $total_count : 0;
		my $freq_N = defined $count_N ? $count_N / $total_count : 0;
		$pos++;
		my $to_print = sprintf ("%d\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.1f\n",$pos,$total_count,$freq_A,$freq_T,$freq_G,$freq_C,$freq_N);
		print $to_print;
}
