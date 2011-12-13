#!/bin/bash
## this is to automate the committing of any new changes to my git repository
echo "----------" >> ~/sync_commit.log
date >> ~/sync_commit.log
echo "---------" >> ~/sync_commit.log
cd ~/src/scripts_by_sofia
rsync -u ~/bin/* ~/src/scripts_by_sofia/.
/usr/local/bin/git add *
/usr/local/bin/git status >> ~/sync_commit.log
/usr/local/bin/git commit -m 'daily commit of any changes to all scritps'
/usr/local/bin/git push
