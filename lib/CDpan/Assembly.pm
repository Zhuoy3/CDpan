#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-04-07

package CDpan::Assembly;

use strict;
use warnings;
use Config::IniFiles;
use File::Copy;

sub assembly {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
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
    #TODO NUM_THREADS是否运行用户指定

    open CONFIG, '>', "$output_dir/$idv_folder_name.masurca_config.txt"
        or die "Error: Couldn't create $output_dir/$idv_folder_name.masurca_config.txt: $!\n";
    print CONFIG $config;
    close CONFIG;

    my $masurca = $par->val('TOOLS', 'masurca');

    my $cmd_masurca = "$masurca $output_dir/$idv_folder_name.masurca_config.txt > $output_dir/$idv_folder_name.masurca_config.log";
    print "Start use cmd: \'$cmd_masurca\'.\n";
    system $cmd_masurca
        and die "Error: Command \'$cmd_masurca\' failed to run normally: $?\n";

    mkdir "$output_dir/assemble"
        or die "Couldn't create $output_dir/assemble: $!\n";
    chdir "$output_dir/assemble";
    system "../assemble.sh > ../assemble.log"
        and die "Error: Command ../assemble.sh: $?\n";
    if ( -e "./CA/final.genome.scf.fasta"){
        move("./CA/final.genome.scf.fasta", "$output_dir/$idv_folder_name.final.genome.scf.fasta")
            or die "Error: Couldn't move file: /CA/final.genome.scf.fasta: $!\n";
    }

    chdir $output_dir;
    # system "rm -rf $output_dir/assemble";
    #TODO暂时关闭
    chdir $main::folder_process;

    return 1;
}


1;

__END__
