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
use CDpan::GetPar;
use CDpan::Check;
use CDpan::QualityControl;
use CDpan::Comparison;
use CDpan::Extract;


our $opt_d = 1;

my $par = new Config::IniFiles(-file => "/storage1/active/zhuoy/test/example.ini");

our $folder_process = "/storage1/active/zhuoy/test/output";

my $idv_folder = "/storage1/active/zhuoy/test/data/sampleA";

our $ref_index = 0; # mark the existence of the reference genome index used by bwa
our $ref_dict = 0; # mark the existence of the reference genome dict used by gatkmy

    my $idv_name = pop [ splitdir($idv_folder) ];
    my $idv_output_folder = catdir($folder_process, $idv_name);

    CDpan::Extract::extract($par, $idv_name, $idv_output_folder)
        or die "ERROR: Operation Extract is abnormal.\n";

