#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpan::Align;

use strict;
use warnings;
use Config::IniFiles;

sub align {
    # &extract($opt, $idv_folder_name, $output_dir, $work_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir, my $work_dir) = @_;

    # Read the software path and set it to the default value
    my $bowtie2 = $par->val('TOOLS', 'bowtie2');
    my $thread = $par->val('BOWTIE2', 'thread');

    my $cmd_align = "$bowtie2 " .
                    "-x $work_dir/index " .
                    "-U $output_dir/$idv_folder_name.singleUnmapped_R2.fq,$output_dir/$idv_folder_name.singleUnmapped_R2.fq " .
                    "-S $output_dir/$idv_folder_name.readContigAlignment.final.sam " .
                    "-p $thread " .
                    "2> $output_dir/$idv_folder_name.align.log";
    print "Start use cmd: \'$cmd_align\'.\n";
    system $cmd_align
        and die "Error: Command \'$cmd_align\' failed to run normally: $?\n";

    return 1;
}

1;

__END__
