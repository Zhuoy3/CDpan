#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-07-21

package CDpan::MMSeqs;

use strict;
use warnings;
use Config::IniFiles;

sub mmseqs {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    # Read the software path and set it to the default value
    my $mmseqs = $par->val('TOOLS', 'mmseqs');

    my $thread = $par->val('MMSEQS', 'thread');

    my $cmd_mmseqs = "$mmseqs easy-linclust " .
                     "$output_dir/$idv_folder_name.filtered.fa " .
                     "$output_dir/$idv_folder_name.filtered.mmseqs " .
                     "$output_dir/mmseqs_tmp " .
                     "--threads $thread " .
                     "--cov-mode 1 " .
                     "-c 0.9 " .
                     "--min-seq-id 0.9 " .
                     "--cluster-mode 2 " .
                     "> $output_dir/$idv_folder_name.filtered.mmseqs.log";
    print "Start use cmd: \'$cmd_mmseqs\'.\n";
    system $cmd_mmseqs
        and die "Error: Command \'$cmd_mmseqs\' failed to run normally: $?\n";

    #TODO 建立索引的程序是否需要

    return 1;
}

1;

__END__
