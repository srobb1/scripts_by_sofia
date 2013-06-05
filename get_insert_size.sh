FQ_1=$1
FQ_2=$2

mkdir tmp-getIS
cd tmp-getIS
head -400000 ../$1 > p1.fq 
head -400000 ../$2 > p2.fq
bwa aln -t 8 -q 20 ~/Wessler-Rice/Genome/index/MSU_r7.all.fa p1.fq > p1.sai 2> out
bwa aln -t 8 -q 20 ~/Wessler-Rice/Genome/index/MSU_r7.all.fa p2.fq > p2.sai 2> out
bwa sampe ~/Wessler-Rice/Genome/index/MSU_r7.all.fa  p1.sai p2.sai p1.fq p2.fq > out.sam 2> out.info 
grep 'inferred external isize from' out.info | perl -pe 's/^.+: //'
cd ..
rm -rf tmp-getIS
