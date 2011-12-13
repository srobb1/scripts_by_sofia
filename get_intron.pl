#!/usr/bin/perl - w
use strict;

my $gff = shift;

my %exons;
my %strand;
my %chr;
open INGFF, $gff or die "can't open $gff $!\n";
while (my $line = <INGFF>){
	chomp $line;
        my @GFF = split /\t/ , $line;
        my ($chr, $type, $start, $end, $strand, $ninth) = ($GFF[0],$GFF[2],$GFF[3],$GFF[4],$GFF[6],$GFF[8]);
	next unless $type eq 'exon';
	#$ninth =~ /Parent=(.+?);/; #if line and id ends with ';'
        $ninth =~ /Parent=(.+?)$/; #if line and id ends with end of line
	my $parent = $1;
        $ninth =~ /ID=(.+?);/; #if line and id ends with ';'
        #$ninth =~ /ID=(.+?)$/; #if line and id ends with end of line
	my $ID = $1;
	$exons{$parent}{$ID}{'start'}=$start;
	$exons{$parent}{$ID}{'end'}=$end;
	$strand{$parent}=$strand;
	$chr{$parent}=$chr;
}

foreach my $parent (sort keys %exons){
	my @coords;
	my $strand =  $strand{$parent};
	my $chr =  $chr{$parent};
	foreach my $exonID (sort {$exons{$parent}{$a}{'start'} <=> $exons{$parent}{$b}{'start'}} keys %{$exons{$parent}}){
		## all exons for a single gene
		my $start =  $exons{$parent}{$exonID}{'start'};
		my $end =  $exons{$parent}{$exonID}{'end'};
		push @coords , ($start , $end);
	}
	#(start, end, start, end, start, end)
	# throw away first and last
	pop @coords;
	shift @coords;
	#now: (end, start, end, start)
	my $intron_total = (scalar @coords) / 2;
	my $intron_count = 1;
	if ($strand eq '-'){
		$intron_count = $intron_total;
	}
	for (my $i=0 ; $i< scalar @coords ; $i=$i+2 ){
		my $intron_start = $coords[$i] + 1;
		my $intron_end = $coords[$i+1] - 1;
		my $intron_ID = $parent.":intron_".$intron_count;
		print "$chr\t.\tintron\t$intron_start\t$intron_end\t.\t$strand\t.\tID=$intron_ID;Parent=$parent\n";
		if ($strand eq '-'){
			$intron_count--;
		}else{
			$intron_count++;
		}		
	}
}
