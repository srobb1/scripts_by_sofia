#!/bin/bash
## this is to automate the committing of any new changes to my git repository
echo "----------" >> sync_commit.log
date >> sync_commit.log
echo "---------" >> sync_commit.log
cd ~/src/scripts_by_sofia
rsync -u ~/bin/* ~/src/scripts_by_sofia/.
git add *
git status >> sync_commit.log
git commit -m 'daily commit of any changes to all scritps'
git push
