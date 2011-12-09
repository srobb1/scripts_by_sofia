#!/usr/bin/perl -w

use strict;
use File::Spec;


##get fq output from sam files

my $dir         = shift;
my $dir_path    = File::Spec->rel2abs($dir);
my $current     = File::Spec->curdir();
my $current_dir = File::Spec->rel2abs($current);

opendir( DIR, $dir ) || die "$!";
foreach my $file ( readdir(DIR) ) {
    my $filename = "$dir_path/$file";
    if ( $file =~ /^(\S+)\.sam$/ ) {
        my $filebase = $1;
        `rm -f $dir_path/$filebase.toMove.sh` if -e "$current_dir/toMove.sh";

        open OUTSH, ">$dir_path/$filebase.sam2fq.sh";
        print OUTSH "#!/bin/bash\n\n";
        print OUTSH "tmp_dir=\`mktemp --tmpdir=/scratch -d\`\n";
        print OUTSH 'cd $tmp_dir', "\n";
        if ( $filebase =~ /unPaired/ ) {
            print OUTSH
"java -jar /home_stajichlab/robb/stajichlab-shared/pkg/Picard/picard/SamToFastq.jar INPUT=$dir_path/$file FASTQ=\$tmp_dir/$filebase.fq VALIDATION_STRINGENCY=LENIENT\n";
        }
        else {
            print OUTSH
"java -jar /home_stajichlab/robb/stajichlab-shared/pkg/Picard/picard/SamToFastq.jar INPUT=$dir_path/$file FASTQ=\$tmp_dir/$filebase"
              . "_1.fq SECOND_END_FASTQ=\$tmp_dir/$filebase"
              . "_2.fq VALIDATION_STRINGENCY=LENIENT\n";
        }
        print OUTSH 'remote_host=`hostname`', "\n";
        print OUTSH "echo \"scp \$remote_host:\$tmp_dir/$filebase"
          . "_1.fq $dir_path/.\" > $current_dir/$filebase.toMove.sh\n";
        print OUTSH "echo \"scp \$remote_host:\$tmp_dir/$filebase"
          . "_2.fq $dir_path/.\" >> $current_dir/$filebase.toMove.sh\n";

        #print OUTSH 'cp $tmp_dir'."/$filebase"."_2.fq $dir_path/.\n";
        print OUTSH "cd $current_dir\n";
        print OUTSH
"echo \"ssh \$remote_host rm -fr \$tmp_dir\" >> $current_dir/$filebase.toMove.sh\n";
        close OUTSH;
    }
}
