#!/usr/bin/perl

##parse primer3 output
##output tab delimited file

use strict;

my $file = shift;
open IN, $file;
my $id;
my $count;
my %primer;
while (my $line = <IN>){
	if ($line =~ /SEQUENCE_ID=(.+)/){
		$id = $1;
		$count = 0;
	}	
	elsif ($line =~ /^SEQUENCE=(.+)/){
		$primer{$id}{1}{SEQUENCE}=$1;
	}
	elsif ($line =~ /PRIMER_PRODUCT_SIZE_RANGE=.+\d+\-(\d+)$/){
		$primer{$id}{1}{SEQ_LENGTH}=$1;
	}
	elsif ($line =~ /PRIMER_LEFT_(\d*_)?SEQUENCE=(.+)/){
		$count++;	
		$primer{$id}{$count}{PRIMER_LEFT_SEQUENCE}=$2;
	}
	elsif ($line =~ /PRIMER_RIGHT_(\d*_)?SEQUENCE=(.+)/){
		$primer{$id}{$count}{PRIMER_RIGHT_SEQUENCE}=$2;;
	}
	elsif ($line =~ /PRIMER_PAIR_PENALTY(_\d*)?=(.+)/){
		$primer{$id}{$count}{PRIMER_PAIR_PENALTY}=$2;;
	}
	elsif ($line =~ /PRIMER_LEFT(_\d*)?=(\d+),(\d+)/){
		$primer{$id}{$count}{PRIMER_LEFT_START}=$2;
		$primer{$id}{$count}{PRIMER_LEFT_LEN}=$3;
	}
	elsif ($line =~ /PRIMER_RIGHT(_\d*)?=(\d+),(\d+)/){
		$primer{$id}{$count}{PRIMER_RIGHT_START}=$2;
		$primer{$id}{$count}{PRIMER_RIGHT_LEN}=$3;
	}
	elsif ($line =~ /PRIMER_LEFT(_\d*)?_TM=(.+)/){
		$primer{$id}{$count}{PRIMER_LEFT_TM}=$2;
	}
	elsif ($line =~ /PRIMER_RIGHT(_\d*)?_TM=(.+)/){
		$primer{$id}{$count}{PRIMER_RIGHT_TM}=$2;
	}
	elsif ($line =~ /PRIMER_LEFT(_\d*)?_GC_PERCENT=(.+)/){
		$primer{$id}{$count}{PRIMER_LEFT_GC}=$2;
	}
	elsif ($line =~ /PRIMER_RIGHT(_\d*)?_GC_PERCENT=(.+)/){
		$primer{$id}{$count}{PRIMER_RIGHT_GC}=$2;
	}
	elsif ($line =~ /PRIMER_PRODUCT_SIZE(_\d*)?=(.+)/){
		$primer{$id}{$count}{PRIMER_PRODUCT_SIZE}=$2;
	}else {
		next;
	}
}
open OUTPRIMERS , ">$file.primersList.txt";
print "ID\tSEQLength\tprimerSetNum\tprimerOrient\tproduct_size\tstart\tlen\ttm\tgc%\tprimerSeq\n";
foreach my $id (sort keys %primer){
	my $seqLen = $primer{$id}{1}{SEQ_LENGTH};
	my $sequence = $primer{$id}{1}{SEQUENCE};		
	foreach my $count (sort {$a <=> $b} keys %{$primer{$id}}){
		next if $count == 0;
                next if $primer{$id}{$count}{PRIMER_LEFT_SEQUENCE} eq '' ; 
		my $leftSeq = $primer{$id}{$count}{PRIMER_LEFT_SEQUENCE};
		my $rightSeq = $primer{$id}{$count}{PRIMER_RIGHT_SEQUENCE};
		my $leftStart = $primer{$id}{$count}{PRIMER_LEFT_START};
		my $leftLen = $primer{$id}{$count}{PRIMER_LEFT_LEN};
		my $rightStart = $primer{$id}{$count}{PRIMER_RIGHT_START};
		my $rightLen = $primer{$id}{$count}{PRIMER_RIGHT_LEN};
		my $leftTM = $primer{$id}{$count}{PRIMER_LEFT_TM};
		my $rightTM = $primer{$id}{$count}{PRIMER_RIGHT_TM}; 
		my $leftGC = $primer{$id}{$count}{PRIMER_LEFT_GC};
		my $rightGC = $primer{$id}{$count}{PRIMER_RIGHT_GC};
		my $primerPairPen = $primer{$id}{$count}{PRIMER_PAIR_PENALTY};
		my $productSize = $primer{$id}{$count}{PRIMER_PRODUCT_SIZE};		

		my $product = substr ($sequence, $leftStart, $productSize);
		print "$id\t$seqLen\t$count\tleft\t$productSize\t$leftStart\t$leftLen\t$leftTM\t$leftGC\t$leftSeq\n";
		print "$id\t$seqLen\t$count\tright\t$productSize\t$rightStart\t$rightLen\t$rightTM\t$rightGC\t$rightSeq\n";
                print OUTPRIMERS "$id-$count\t$leftSeq\t$rightSeq\n";
		#print OUTFASTA ">$id(primerSet:$count) productSize=$productSize leftTM=$leftTM leftGC=$leftGC leftSeq=$leftSeq rightTM=$rightTM rightGC=$rightGC rightSeq=$rightSeq\n$product\n";
	}
}
