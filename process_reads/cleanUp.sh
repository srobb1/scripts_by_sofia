if [ !-d trash ] ; then
  mkdir trash
fi
 
for i in `ls | grep -v trash | grep -v flagstat | grep -v recal.ba | grep -v vcf | grep -v tab | grep -v clean` ; do mv $i trash ; done

