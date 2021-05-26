#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-15

package CDpan::Comparison;

use strict;
use warnings;
use Config::IniFiles;
use File::Spec::Functions qw /:ALL/;

sub comparison {
    # &qualitycontrol($opt, $idv_folder_name, $output_dir, $idv_merge_file)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder is a directory of individual containing genomic data
    # $idv_merge_file is a quotation of a list of two files which is merge data of genomic data
    # $idv_merge_file is returned by CDpan::QualityControl::QualityControl

    (my $par, my $idv_folder_name, my $output_dir, my $idv_merge_file) = @_;
    (my $idv_merge_file1, my $idv_merge_file2) = @$idv_merge_file;

    my $output = catfile($output_dir, "${idv_folder_name}.sam" );
    print "bwa output file: $output.\n" if $main::opt_d;

    # Reference genome file is checked in CDpan::Check

    my $thread = $par->val('COMPARISON', 'thread');
    my $ref    = $par->val('COMPARISON', 'ref');

    # Generate reference genome index
    unless ($main::bwa_ref_index) {
        my $cmd_bwa_index = "bwa index $ref";
        print "Start use cmd: \'$cmd_bwa_index\'\n.";
        system $cmd_bwa_index;
        $main::bwa_ref_index = 1;
    }

    #TODO 考虑将所有软件集中至统一模块进行调整
    # Read the software path and set it to the default value
    (my $bwa = $par->val('TOOLS', 'bwa', './bwa') ) =~ s/\/bwa$//;
    # Add environment variables
    $ENV{PATH} = "$bwa:$ENV{PATH}:";

    my $cmd_bwa = "bwa mem -t $thread -M " .
                  "-R \"\@RG\\tID:${idv_folder_name}\\tLB:${idv_folder_name}\\tPL:ILLUMINA\\tSM:${idv_folder_name}\" " .
                  "$ref $idv_merge_file1 $idv_merge_file2 > $output";

    print "Start use cmd: \'$cmd_bwa\'\n.";
    system $cmd_bwa;

    return $output;
}


1;

__END__
