#!/usr/bin/perl -w
use strict;

my $insert_files = shift;
my $existing_files=shift;
my @insert_files;
if (defined $insert_files){
  my @insert_files = split ',' , $insert_files;
}
my @existing_files;
if (defined $exisiting_files){
  my @existing_files = split ',' , $existing_files;
}

