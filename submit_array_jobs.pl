#!/usr/bin/perl -w

use strict;

my $file_of_array_jobs = shift;
open JOBS, "$file_of_array_jobs" or die "Can't open $file_of_array_jobs";

while (my $job = <JOBS>){
  chomp $job;
  die "not a qsub array job, $job" if $job !~ /qsub\s+\-t/;
  print "checking to see if i can submit: $job\n";
  my $empty_Q = 0;
  while (!$empty_Q){
    my $status = `qstat | grep robb`;
    chomp $status;
    print "status: $status\n";
    if ($status){
      `sleep 5m`;
    }else{
      $empty_Q = 1;
      print "submitting next job: $job\n";
      `$job`;      
    }
  }
}
