#!/usr/bin/perl -w
use strict;


## lengths_file OGs_file @mod_files

###==> concat_OGs.txt <==
###>Allo
###(OG1035)(OG1044)(OG1083)(OG1085)(OG1114)(OG1212)(OG1223)(OG1246)(OG1270)(OG1343)(OG1383)(OG1388)(OG1422)(OG1425)(OG1487)(OG1630)(OG1634)(OG1700)(OG1827)(OG1920)(OG1935)(OG2000)(OG2015)(OG2016)(OG2091)(OG2172)(OG2255)(OG2278)(OG2295)(OG2327)(OG2660)(OG2921)(OG2992)(OG3080)(OG3090)(OG3507)(OG3516)(OG3751)(OG406)(OG665)(OG686)(OG879)(OG994)
###
###==> concat_lens.txt <==
###>Allo
###(827)(461)(306)(622)(479)(542)(592)(476)(242)(348)(462)(489)(540)(192)(777)(390)(441)(81)(210)(250)(431)(433)(383)(410)(348)(209)(236)(257)(294)(788)(210)(231)(157)(296)(254)(264)(340)(181)(446)(233)(1317)(560)(359)

###==> ModelTest/OG1035.pep.19.trimal.aln.phy.modelselection.out <==
###Best Model : RTREVF


##to print this
#DNA, p1=1-30
#DNA, p2=31-60

my $lengths_file =shift;
my $OGs_file = shift;
my @mod_files = @ARGV;
my %mods;

foreach my $mod_file(@mod_files){
  my $line = `head -n 1 $mod_file`;
  my ($OG) = $mod_file =~ /(OG\w+)\./;
  my ($mod) = $line =~ /Best.Model.:.(.+)$/;
  $mods{$OG}=$mod;
}


my $lengths_line = `head -n2 $lengths_file | tail -n1`;
my $OGs_line = `head -n2 $OGs_file | tail -n1`;
my (@lens) = $lengths_line =~ /\((\d+)\)/g;
my (@OGs) = $OGs_line =~ /\((\w+)\)/g;

my $total = 0;
my $i = 1;
for (my $j = 0; $j < @lens ; $j++){
  my $len = $lens[$j];
  my $OG = $OGs[$j];
  my $mod = $mods{$OG};
  my $start = $total +1;
  $total = $total + $len;
  print "$mod, p$i=$start-$total\n";
  $i++; 
}
