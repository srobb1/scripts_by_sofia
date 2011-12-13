#!/usr/bin/perl -w
use strict;
use URI::Escape;
use Data::Dumper;

# ./get_intron.pl MSU_r7.all.gff3 > introns.gff
# ./get_intergenic.pl MSU_r7.all.gff3 > intergenic.gff
# cat MSU_r7.all.gff3 introns.gff intergenic.gff > all_plus_intergenic_introns.gff 
# perl -pi -e 's/^\s+$//' all_plus_intergenic_introns.gff
# ./combine_data.pl Chr1.mPing.te_insertion_sites.table.txt all_plus_intergenic_introns.gff all.GOSlim_assignment > out.txt

my $insert = shift; #table ouput file from find_mping_insertions.pl
my $gff	   = shift; #from MSU of all genes, exons, and UTR
my $go	   = shift; #from MSU of all GOterms for genes
my %inserts;
my %annot;
my %go;
my %tds;
open (INTABLE, "<", $insert) or die "Can't open $insert $!\n";
<INTABLE>;
while (my $line = <INTABLE>){
	chomp $line;
#	my ($te, $chr, $start, $end, $four, $left_flanking_seq, ) = split /\t/ , $line;

	my ($chr, $start, $end, $four, $left_flanking_seq, ) = split /\t/ , $line;
	$tds{$chr}{$start}= substr($left_flanking_seq, -3);
	push @{$inserts{$chr}} , $start;
}
close INTABLE;

open (GO, "<" , $go) or die "Can't open $go $! \n";
<GO>;
while (my $line = <GO>){
	chomp $line;
	my ($loc, $goid, $goterm, $gopart)= split /\t/, $line;
	$loc =~ s/.\d+$//;
	$go{$loc}{$goid}{'part'}=$gopart;
	$go{$loc}{$goid}{'term'}=$goterm;
}
close GO;

open (INGFF, "<", $gff) or die "Can't open $gff $!\n";
<INGFF>;
while (my $line = <INGFF>){
	chomp $line;
	my @GFF = split /\t/ , $line;
	my ($chr, $type, $start, $end, $strand, $nine) = ($GFF[0], $GFF[2] , $GFF[3], $GFF[4], $GFF[6], $GFF[8]);
	next if !exists $inserts{$chr};
	next if $type eq 'CDS' or $type eq 'mRNA';
	$nine =~ /ID=(.+?);/;
	my $ID = $1;
	
	foreach my $insert (sort {$a <=> $b} @{$inserts{$chr}}){
		if ($start <= $insert and $insert <= $end){
		#insert is in this feature
		$annot{$chr}{"$start..$end"}{'start'}=$start;
		$annot{$chr}{"$start..$end"}{'end'}=$end;
		$annot{$chr}{"$start..$end"}{'type'}=$type;
		$annot{$chr}{"$start..$end"}{'ID'}=$ID;
		$annot{$chr}{"$start..$end"}{'insert'}=$insert;
		$annot{$chr}{"$start..$end"}{'strand'}=$strand;
			if ($type eq 'gene'){
				$nine =~ /Note=(.+)$/;
				my $desc = $1;
				my $str  = uri_unescape($desc);
				$annot{$chr}{"$start..$end"}{'desc'}=$str;
			}elsif ($type eq 'intergenic'){
				$nine =~ /downstream_gene=(.+?);/;
				$annot{$chr}{"$start..$end"}{'downstream'} = $1;
				$nine =~ /upstream_gene=(.+?)$/;
				$annot{$chr}{"$start..$end"}{'upstream'} = $1;
				
			}
		}		
	}
}
close INGFF;


print "chr\tinsert_pos\tref TDS\tfeature\tstrand\tf.start\tf.end\tf.ID\tf.Desc\tintergenic.down.gene\tintergenic.up.gene\tgene.go.id\tgene.go.part\tgene.go.term\n";		
foreach my $chr  (keys %annot){
	foreach my $range ( sort { (split /\.\./, $a)[0] <=> (split /\.\./, $b)[0]} keys %{$annot{$chr}}){
		my $start = $annot{$chr}{$range}{'start'};
		my $end = $annot{$chr}{$range}{'end'};
		my $type = $annot{$chr}{$range}{'type'};
		my $insert = $annot{$chr}{$range}{'insert'};
		my $ID = $annot{$chr}{$range}{'ID'};
		my $strand = $annot{$chr}{$range}{'strand'};
		if ($strand eq '-'){
			($start , $end) = ($end , $start);
		}
		my $note = 'N/A';
		my $downstream = 'N/A';
		my $upstream = 'N/A';
		my $go = '';
		if ($type eq 'gene'){
			$note = $annot{$chr}{$range}{'desc'};
			if (exists $go{$ID}){
				foreach my $go_id (keys %{$go{$ID}}){ 
					my $goterm = $go{$ID}{$go_id}{'term'};
					my $part = $go{$ID}{$go_id}{'part'};
					$go .= "$go_id\t$part\t$goterm\t";
				}
			}		
		}elsif($type eq 'intergenic'){
			$downstream = $annot{$chr}{$range}{'downstream'};
			$upstream = $annot{$chr}{$range}{'upstream'};
			$ID = 'N/A';
		}
		my $tds = $tds{$chr}{$insert};

		print "$chr\t$insert\t$tds\t$type\t$strand\t$start\t$end\t$ID\t$note\t$downstream\t$upstream\t$go\n";		
	}
}
