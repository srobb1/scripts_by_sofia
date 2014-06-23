cd /rhome/robb/Wessler-Rice/RIL/genotypes/MSU_r7.corrected

if [ ! -d trash ] ; then 
  mkdir trash
fi

  mv *.realign.bai trash
  mv *.realign.bam trash
  mv *.dedup.metrics trash
  mv *.dedup.bai trash
  mv *.dedup.bam trash
  mv *.recal_data.grp trash
  mv *.sam trash
  mv *.intervals trash
  mv *_p2.sai trash
  mv *_p1.sai trash
  mv *.RG.bai trash
  mv *.RG.bam trash
