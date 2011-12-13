#!/usr/bin/perl - w
use strict;

my $gff = shift;
##edit to be specific for your organism
my %chr = ('Chr1',43270923,'Chr2',35937250,'Chr3',36413819,'Chr4',35502694,'Chr5',29958434,'Chr6',31248787,'Chr7',29697621,'Chr8',28443022,'Chr9',23012720,'Chr10',23207287,'Chr11',29021106,'Chr12',27531856);

open INGFF, $gff or die "can't open $gff $!\n";
my ($upstream ,$last_upstream ,$downstream) = ('','N/A','');
my $last_end = 0;
#change this to correspond to your first Chr in the gff
my $last_chr='Chr1';
while (my $line = <INGFF>){
        chomp $line;
        my @GFF = split /\t/ , $line;
        my ($chr, $type, $start, $end, $strand, $ninth) = ($GFF[0],$GFF[2],$GFF[3],$GFF[4], $GFF[6],$GFF[8]);
        next unless $type eq 'gene';
	$ninth =~ /ID=(.+?);/; #if line and id ends with ';'
	#$ninth =~ /ID=(.+?)$/; #if line and id ends with end of line
	$upstream = "$1:$start..$end($strand)";
	$downstream = $last_upstream;
        my ($intergenic_start,$intergenic_end) = ( $last_end+1 , $start-1 );
        if ($last_chr ne $chr){
                $last_end = 0;
                my $chr_end = $chr{$last_chr};
		$last_upstream = 'N/A';
        	print "$last_chr\t.\tintergenic\t$intergenic_start\t$chr_end\t.\t.\t.\tID=$chr.intergenic.$intergenic_start.$chr_end;downstream_gene=$downstream;upstream_gene=N/A\n" if $intergenic_start < $chr_end;
		$intergenic_start = 1;
                $last_chr=$chr;
        }
        print "$chr\t.\tintergenic\t$intergenic_start\t$intergenic_end\t.\t.\t.\tID=$chr.intergenic.$intergenic_start.$intergenic_end;downstream_gene=$downstream;upstream_gene=$upstream\n" if $intergenic_start < $intergenic_end;
        $last_end = $end;
	$last_upstream = $upstream;
}
