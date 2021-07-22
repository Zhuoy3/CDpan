#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-25

package CDpan::Extract;

use strict;
use warnings;
use Config::IniFiles;

sub extract {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output

    (my $par, my $idv_folder_name, my $output_dir) = @_;

    # Read the software path and set it to the default value
    my $samtools = $par->val('TOOLS', 'samtools');

    my $thread = $par->val('EXTRACT', 'thread');

    my $cmd_fastq1 = "$samtools fastq -\@ $thread " .
                     "-f 12 $output_dir/$idv_folder_name.sort.bam " .
                     "-1 $output_dir/$idv_folder_name.mateUnmapped_R1.fq " .
                     "-2 $output_dir/$idv_folder_name.mateUnmapped_R2.fq " .
                     "2> /dev/null";
    print "Start use cmd: \'$cmd_fastq1\'.\n";
    system $cmd_fastq1
        and die "Error: Command \'$cmd_fastq1\' failed to run normally: $?.\n";

    my $cmd_fastq2 = "$samtools fastq -\@ $thread " .
                     "-f 68 -F 8 $output_dir/$idv_folder_name.sort.bam " .
                     "> $output_dir/$idv_folder_name.R1_mateMapped.fq " .
                     "2> /dev/null";
    print "Start use cmd: \'$cmd_fastq2\'.\n";
    system $cmd_fastq2
        and die "Error: Command \'$cmd_fastq2\' failed to run normally: $?.\n";

    my $cmd_fastq3 = "$samtools fastq -\@ $thread " .
                     "-f 132 -F 8 $output_dir/$idv_folder_name.sort.bam " .
                     "> $output_dir/$idv_folder_name.R2_mateMapped.fq " .
                     "2> /dev/null";
    print "Start use cmd: \'$cmd_fastq3\'.\n";
    system $cmd_fastq3
        and die "Error: Command \'$cmd_fastq3\' failed to run normally: $?.\n";

    my $cmd_view = "$samtools view -\@ $thread " .
                   "-f 8 -F 4 -b -h $output_dir/$idv_folder_name.sort.bam " .
                   "> $output_dir/$idv_folder_name.Sus_11_1_Links.bam ";
    print "Start use cmd: \'$cmd_view\'.\n";
    system $cmd_view
        and die "Error: Command \'$cmd_view\' failed to run normally: $?.\n";

    my $cmd_fastq4 = "$samtools fastq -\@ $thread " .
                     "-f 4 $output_dir/$idv_folder_name.sort.bam " .
                     "-1 $output_dir/$idv_folder_name.singleUnmapped_R1.fq " .
                     "-2 $output_dir/$idv_folder_name.singleUnmapped_R2.fq " .
                     "2> /dev/null";
    print "Start use cmd: \'$cmd_fastq4\'.\n";
    system $cmd_fastq4
        and die "Error: Command \'$cmd_fastq4\' failed to run normally: $?.\n";

    return 1;
}


1;

__END__
