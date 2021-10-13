#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-07-21

package CDpan::Nucmer;

use strict;
use warnings;
use Config::IniFiles;

sub nucmer {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    # Read the software path and set it to the default value
    my $nucmer = $par->val('TOOLS', 'nucmer');

    my $qry = $par->val('DATA', 'qry');

    my $cmd_nucmer = "$nucmer " .
                     "-p $output_dir/$idv_folder_name.filtered.mmseqs " .
                     "$qry " .
                     "$output_dir/$idv_folder_name.filtered.mmseqs_rep_seq.fasta " .
                     "-g 1000 " .
                     "-c 90 " .
                     "-l 500 " .
                     "-t 10";
    print "Start use cmd: \'$cmd_nucmer\'.\n";
    system $cmd_nucmer
        and PrintErrorMessage("Command \'$cmd_nucmer\' failed to run normally: $?\n");

    my $show_coords = $par->val('TOOLS', 'show-coords');

    my $cmd_show_coords = "$show_coords " .
                          "-rcl " .
                          "$output_dir/$idv_folder_name.filtered.mmseqs.delta " .
                          "> $output_dir/$idv_folder_name.filtered.mmseqs.delta.coords";
    print "Start use cmd: \'$cmd_show_coords\'.\n";
    system $cmd_show_coords
        and PrintErrorMessage("Command \'$cmd_show_coords\' failed to run normally: $?\n");

    return 1;
}

1;

__END__
