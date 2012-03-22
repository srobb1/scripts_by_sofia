##shell script to check the log file of raw_reads_process_clean
## ex: check_raw_log.sh 2012-Mar-21-log.txt

log=$1

echo "log parsed:"
awk {'print $1'} $log | sort | uniq -c | sort -n
echo ""
echo "log line count:"
awk {'print $1'} $log | sort | uniq -c | sort -n | wc -l 
echo ""
echo "2 count => good:"
awk {'print $1'} $log | sort | uniq -c | awk {'print $1'} | grep -c 2
echo "1 count => bad, you want 0 1 counts:"
awk {'print $1'} $log | sort | uniq -c | awk {'print $1'} | grep -c 1 

