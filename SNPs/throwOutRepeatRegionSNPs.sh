vcf=$1
genome_repeat=$2
base=`basename $vcf .raw.vcf`
cwd=`pwd`
grep ^# $vcf > file_noRepeats.vcf
subtractBed -a $vcf -b $genome_repeat >> file_noRepeats.vcf
perl -pi -e 's/\t$//' file_noRepeats.vcf
mv file_noRepeats.vcf $cwd/$base.noRepeats.vcf
