#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-08-21

package CDpan::Integration;

use strict;
use warnings;
use Config::IniFiles;

sub integration {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    mkdir $output_dir/$link_new or die "Error: Cannot create process folder '$output_dir/$link_new': $!\n";

    return 1;
}

1;

__END__
