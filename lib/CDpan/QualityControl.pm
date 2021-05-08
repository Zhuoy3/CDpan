#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-02

package CDpan::QualityControl;

use strict;
use warnings;
use File::Spec::Functions;
use Config::IniFiles;

sub qualitycontrol {
    # &qualitycontrol($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    my $par = shift;

    # Get individual folders that need to be processed
    # sample A, B, C, D
        # Get individual files that need to be processed
        # 1, 2, 3
            # processed
    system "";
}
