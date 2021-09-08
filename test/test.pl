#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-03

use strict;
use warnings;
use feature qw /say/;

use FindBin;
use lib "$FindBin::Bin/lib";

use File::Spec::Functions  qw /:ALL/;
use File::Slurp;
use Config::IniFiles;
use CDpan::RepeatMasker;

my $par = new Config::IniFiles(-file => "/storage1/active/zhuoy/CDpan/example.ini");

my $folder_process = "/storage1/active/zhuoy/CDpan/tmp/";

CDpan::RepeatMasker::repeat_masker($par, $folder_process)
    or die "Error: Operation RepeatMasker is abnormal.\n";
