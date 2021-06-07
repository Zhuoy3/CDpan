#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-04-07

package CDpan::Assembly;

use strict;
use warnings;
use Config::IniFiles;

sub assembly {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a directory path which is used to output
    #TODO 未调试
    (my $par, my $idv_folder_name, my $output_dir) = @_;
    chdir $output_dir;

    my $config = "" .
    "DATA\n" .
    "PE= pe 300 50 $output_dir/$idv_folder_name.mateUnmapped_R1.fq $output_dir/$idv_folder_name.mateUnmapped_R2.fq\n" .
    "PE= s1 300 50 $output_dir/$idv_folder_name.R1_mateMapped.fq\n" .
    "PE= s2 300 50 $output_dir/$idv_folder_name.R2_mateMapped.fq\n" .
    "END\n" .
    "\n" .
    "PARAMETERS\n" .
    "GRAPH_KMER_SIZE = auto\n" .
    "USE_LINKING_MATES = 1\n" .
    "KMER_COUNT_THRESHOLD = 1\n" .
    "NUM_THREADS = 24\n" .
    "JF_SIZE=200000000\n" .
    "DO_HOMOPOLYMER_TRIM=0\n" .
    "END\n";

    open CONFIG, '>', "$output_dir/$idv_folder_name.config.txt";
    print CONFIG $config;
    close CONFIG;

    (my $masurca = $par->val('TOOLS', 'masurca', './masurca') ) =~ s/\/masurca$//;
    $ENV{PATH} = "$masurca:$ENV{PATH}:";

    system "masurca $output_dir/$idv_folder_name.config.txt";
    system "";

    chdir $main::folder_process;

    return 1;
}


1;

__END__
