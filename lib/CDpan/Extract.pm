#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-17

package CDpan::Extract;

use strict;
use warnings;
use Config::IniFiles;

sub extract {
    # &qualitycontrol($opt, $idv_folder, $idv_merge_file)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder is a directory of individual containing genomic data
    # $idv_merge_file is a quotation of a list of two files which is merge data of genomic data
    # $idv_merge_file is returned by CDpan::QualityControl::QualityControl

    (my $par, my $idv_folder, my $idv_merge_file) = @_;
    (my $idv_merge_file1, my $idv_merge_file2) = @$idv_merge_file;



    return $output;
}


1;

__END__
