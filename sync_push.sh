#!/bin/bash
## this is to automate the committing of any new changes to my git repository

cd ~/src/scripts_by_sofia
rsync -u ~/bin/* ~/src/scripts_by_sofia/.
git add *
git status >> sync_push.log
git commit -m 'daily commit of any changes to all scritps'
git push
