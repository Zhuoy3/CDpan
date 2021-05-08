#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-02

package CDpan::QualityControl;

use strict;
use warnings;
use File::Spec::Functions qw /:ALL/;
use Config::IniFiles;

sub QualityControl {
    # &qualitycontrol($opt, $idv_folder, $idv_file)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder is a directory of individual containing genomic data
    # $idv_file is a quotation of a list of some files which is genomic data
    (my $par, my $idv_folder, my $idv_file) = @_;
    my @idv_file = @$idv_file;

    # The default value is set when the parameter is imported, and there is no need to handle it here
    my $quality    = $par->val('QUALITYCONTROL', 'quality'); # Trim low-quality ends from reads in addition to adapter removal.
    my $length     = $par->val('QUALITYCONTROL', 'length'); # Discard reads that became shorter than length INT because of either quality or adapter trimming.
    my $error_rate = $par->val('QUALITYCONTROL', 'error_rate'); # Maximum allowed error rate (no. of errors divided by the length of the matching region)
    my $cores      = $par->val('QUALITYCONTROL', 'cores'); # Number of cores to be used for trimming

    my @idv_folder_split = splitdir($idv_folder);
    my $idv_folder_name = pop @idv_folder_split;
    my $output = catdir($main::folder_process, $idv_folder_name);

    # Read the software path and set it to the default value
    (my $trim_galore = $par->val('TOOLS', 'trim_galore', './trim_galore') ) =~ s/\/trim_galore$//;
    (my $cutadapt    = $par->val('TOOLS', 'cutadapt',    './cutadapt')    ) =~ s/\/cutadapt$//;
    (my $fastqc      = $par->val('TOOLS', 'fastqc',      './fastqc')      ) =~ s/\/fastqc$//;

    # Add environment variables
    $ENV{PATH} = "$ENV{PATH}:$trim_galore:$cutadapt:$fastqc:";

    while (@idv_file) {
        my $paired1 = shift @idv_file;
        my $paired2 = shift @idv_file;
        # trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired a1_1.fq.gz a1_2.fq.gz --gzip -o ./ -j 6
        my $cmd_trim_galore = "trim_galore " .
                           "--quality $quality " .
                           "--phred33 " .
                           "--stringency 3 " .
                           "--length $length " .
                           "-e $error_rate " .
                           "--paired $paired1 $paired2 " .
                           "--gzip " .
                           "--output_dir $output " .
                           "-j $cores";

        print "Start use cmd: \'$cmd_trim_galore\'.";
        system 'echo $PATH';
        system $cmd_trim_galore;
    }

    return 1;
}


1;

__END__
