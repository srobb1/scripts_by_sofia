ls -l *.sh.e* | awk {'print $5'} | sort | uniq -c
