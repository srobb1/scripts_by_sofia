#!/bin/bash
## this is to automate the committing of any new changes to my git repository
echo "----------" >> ~/sync_commit.log
date >> ~/sync_commit.log
echo "---------" >> ~/sync_commit.log
cd ~/src/scripts_by_sofia

##scripts_by_sofia
rsync -uL ~/bin/* ~/src/scripts_by_sofia/.
rsync -uL ~/bin/findPackMule/* ~/src/scripts_by_sofia/findPackMule/.
rsync -uL ~/bin/find_mping_insertions/* ~/src/scripts_by_sofia/find_mping_insertions/.
rsync -uL ~/bin/gff_and_seqfeaturestore_scripts/* ~/src/scripts_by_sofia/gff_and_seqfeaturestore_scripts/.
rsync -uL ~/bin/process_raw_paired_reads_2_split_by_target/* ~/src/scripts_by_sofia/process_raw_paired_reads_2_split_by_target/.
rsync -uL ~/bin/dynGenomClass/* ~/src/scripts_by_sofia/dynGenomClass/.

##TEamRice
rsync -uL ~/bin/find_mping_insertions/* ~/src/TEamRice/find_mping_insertions/.
rsync -uL ~/bin/gff_and_seqfeaturestore_scripts/* ~/src/TEamRice/gff_and_seqfeaturestore_scripts/.
rsync -uL ~/bin/process_raw_paired_reads_2_split_by_target/* ~/src/TEamRice/process_raw_paired_reads_2_split_by_target/.


/usr/local/bin/git add *
/usr/local/bin/git status >> ~/sync_commit.log
/usr/local/bin/git commit -m 'daily commit of any changes to all scritps'
/usr/local/bin/git push
