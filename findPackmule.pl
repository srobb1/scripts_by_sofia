#!/usr/bin/perl -w
use strict;
use Getopt::Long;

my $mule_fasta = shift or &getHelp;
my $db_file = shift or &getHelp; #same file used as a database in the above blast
my $db2 = shift or &getHelp; # nr, species specific of the above organsims

my $min_insert = 2000;
my $max_insert = 8000;
my $min_ident  = 85;
my $query_half = 0.4;
my $dir = '~/bin/findPackMule';

GetOptions(
    'h|help'  => \&getHelp,
    'm|min_insert:i' => \$min_insert,
    'x|max_insert:i' => \$max_insert,
    'i|min_ident:i'  => \$min_ident,
    'p|query_half:f' => \$query_half,
);

sub getHelp {
    print "

-------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------
usage:
./find_packmule.pl mule.nt.fa genome.fa species.nt.fa -m 2000 -x 8000 -i 85 -p 0.4 [-h] 

-------------------------------------------------------------------------------------------
options:
-------------------------------------------------------------------------------------------
-h this message
-m INT min insert size (default 2000)
-x INT max insert size (default 8000)
-i INT min blast %identity (85)
-p FRAC the proportion of the TSD-TIR-TIR-TSD that must be found in a single hsp (0.4)


-------------------------------------------------------------------------------------------
mule.nt.fa 
-------------------------------------------------------------------------------------------
a fasta file that needs to contain the TSD-TIR1-TIR2-TSD of the mule elements


sample fasta  
>mule_autonomous TSD=CTTCAAATG
CTTCAAATGGGGTCTACCCCGTTTGGCATAATGCCGTTTGGCATAATGCCGTTTGGCATACAGTCGTTTGGCATAAAGTCGTTTGGC
ATAATAGTCATTTGGCATAACAGTCGTTTGGCATAATGGTCATTTGGCATAATGGTCGTTTGGCATAATTATGCCAAACGACTATTA
TGCCAAATGACCATTATGCCAAATGACTATTATGCCAAATGGCATTATGCCAAACGACTATTATGCCAAACGACTGTATGCCAAACG
GCATTATGCCAAACGGCATTATGCCAAACGGGGTAGACCCCTTCAAATG

-------------------------------------------------------------------------------------------
genome.fa
-------------------------------------------------------------------------------------------
a fasta file of the genome or supercontigs of the species you want to search in
this file can be blast formatted, but if not, the script will format it 

-------------------------------------------------------------------------------------------
species.nt.fa
-------------------------------------------------------------------------------------------

a fasta of nt sequences to use as a blast database to search the packmule inserts against
this file can be blast formatted, but if not, the script will format it 

";
    exit 1;
}


my $blast = "mule_vs_genome.blastout";

if (!-e "$db_file.nin"){
  print "formating $db_file for blast\n\n";
  `formatdb -i $db_file -p F -o T`;
}

if (!-e $blast or -z $blast){
  print "Running blastn of $mule_fasta against $db_file\n\n";
  `blastall -p blastn -d $db_file -i $mule_fasta -o $blast`;
}else{
  print "Found: $blast already exists.  Will use this file to search\n";
}
print "Finding potential PackMule elements, this might take approx 5 minutes or more\n\n";
`perl $dir/parseBlast_mule.pl $blast $db_file $db2 $min_insert $max_insert $min_ident $query_half > parseBlast_mule.generalOutput.txt`;

if (!-e "$db2.nin"){
  print "formating $db_file for blast\n\n";
  `formatdb -i $db_file -p F -o T`;
}
my $inserts_blastout = "inserts_vs_species.blastout";
if (!-e $inserts_blastout or -z $inserts_blastout){ 
  print "Running blastn of inserts against $db2\n\n";
  `blastall -p blastn -d $db2 -i insertOnly.fa -o $inserts_blastout`;
}else{
  print "Found: $inserts_blastout already exists.  Will use this file to search\n";
}

print "Parsing blast results\n\n";
`perl $dir/parseBlast_packmule_inserts.pl $inserts_blastout > inserts_blastHits.table.txt`;


print "Output files:
blast output of mules vs genome:                   mule_vs_genome.blastout
table of indiviual TSD-TIR hits:                   match.table.txt
fasta containing inserts:                          insertsOnly.fa
fasta containing potential packmule elements:      packmule.fa
blast output of inserts vs species nt seqs:        inserts_vs_species.blastout
table of insert blast hits:                        inserts_blastHits.table.txt  	
";
