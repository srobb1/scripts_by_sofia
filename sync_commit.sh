#!/bin/bash
## this is to automate the committing of any new changes to my git repository

echo "----------" >> ~/sync_commit.log
date >> ~/sync_commit.log
echo "---------" >> ~/sync_commit.log

##scripts_by_sofia


rsync -uL ~/bin/* ~/src/scripts_by_sofia/.
rsync -uL ~/bin/findPackMule/* ~/src/scripts_by_sofia/findPackMule/.
rsync -uL ~/bin/find_mping_insertions/* ~/src/scripts_by_sofia/find_mping_insertions/.
rsync -uL ~/bin/gff_and_seqfeaturestore_scripts/* ~/src/scripts_by_sofia/gff_and_seqfeaturestore_scripts/.
rsync -uL ~/bin/process_raw_paired_reads_2_split_by_target/* ~/src/scripts_by_sofia/process_raw_paired_reads_2_split_by_target/.
rsync -uL ~/bin/dynGenomClass/* ~/src/scripts_by_sofia/dynGenomClass/.
rsync -uL ~/bin/relocaTE/* ~/src/scripts_by_sofia/relocaTE/.
rsync -uL ~/bin/relocaTE_editing/* ~/src/scripts_by_sofia/relocaTE_editing/.
rsync -uL ~/bin/maize/* ~/src/scripts_by_sofia/maize/.
rsync -uL ~/bin/SNPs/* ~/src/scripts_by_sofia/SNPs/.
rsync -uL ~/bin/pong/* ~/src/scripts_by_sofia/pong/.
rsync -uL ~/bin/andrii/* ~/src/scripts_by_sofia/andrii/.
rsync -uL ~/bin/emptySites/* ~/src/scripts_by_sofia/emptySites/.

##TEamRice
#rsync -uL ~/bin/find_mping_insertions/* ~/src/TEamRice/find_mping_insertions/.
#rsync -uL ~/bin/gff_and_seqfeaturestore_scripts/* ~/src/TEamRice/gff_and_seqfeaturestore_scripts/.
#rsync -uL ~/bin/process_raw_paired_reads_2_split_by_target/* ~/src/TEamRice/process_raw_paired_reads_2_split_by_target/.


cd ~/src/scripts_by_sofia
/usr/bin/git add *
/usr/bin/git status >> ~/sync_commit.log
/usr/bin/git commit -m 'daily commit of any changes to all scritps'
/usr/bin/git push

#cd ~/src/TEamRice
#/usr/bin/git add *
#/usr/bin/git status >> ~/sync_commit.log
#/usr/bin/git commit -m 'daily commit of any changes to all scritps'
#/usr/bin/git push
