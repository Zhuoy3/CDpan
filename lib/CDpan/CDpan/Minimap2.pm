#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-08-21

package CDpan::Integration;

use strict;
use warnings;
use Config::IniFiles;
use File::Copy qw / copy move /;
use List::Util qw / max min /;
use File::Path qw / mkpath /;

sub integration {
    # &extract($opt, $idv_folder_name, $output_dir, $work_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir, my $work_dir) = @_;

    my $minimap2 = $par->val('TOOLS', 'minimap2');

    mkdir "$output_dir/minimap2"
        or PrintErrorMessage("Cannot create process folder '$output_dir/minimap2': $!\n");
    chdir "$output_dir/minimap2"
        or PrintErrorMessage("Cannot chdir to '$output_dir/minimap2': $!\n");

    my $cmd_minimap2 = "$minimap2 " .
                       "-x asm10 " .
                       "$work_dir/all.fasta.fai " .
                       "$output_dir/$idv_folder_name.filtered.mmseqs.final.fa " .
                       "> $output_dir/minimap2/aln.paf";
    print "Start use cmd: \'$cmd_minimap2\'.\n";
    system $cmd_minimap2
        and PrintErrorMessage("Command \'$cmd_minimap2\' failed to run normally: $?\n");

    return 1;
}


1;

__END__
