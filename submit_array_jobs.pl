#!/usr/bin/perl -w

use strict;

my $file_of_array_jobs = shift;
open JOBS, "$file_of_array_jobs" or die "Can't open $file_of_array_jobs";
my $usrname = `whoami`;
chomp $usrname;
while (my $job = <JOBS>){
  chomp $job;
  die "not a qsub job, $job" if $job !~ /qsub/;
  #warn "not a qsub array job, $job" if $job !~ /qsub\s+\-t/;
  print "checking to see if i can submit: $job\n";
  my $empty_Q = 0;
  while (!$empty_Q){
    my $status;
    my @status = `qstat | grep $usrname`;
    chomp @status;
    if (@status > 1){
      $status = join "\n" , @status;
    }elsif (@status) {
      $status = $status[0];
    }
    print "status: $status\n" if defined $status;
    if ($status){
      `sleep 5m`;
    }else{
      $empty_Q = 1;
      print "submitting next job: $job\n";
      `$job`;      
    }
  }
}
