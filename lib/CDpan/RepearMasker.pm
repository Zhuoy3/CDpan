#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-08-11

package CDpan::RepeatMasker;

use strict;
use warnings;
use Config::IniFiles;

sub repeat_masker {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    # Read the software path and set it to the default value
    my $repeat_masker = $par->val('TOOLS', 'RepeatMasker');

    my $cmd_repeat_masker = "$repeat_masker " .
                            "-nolow " .
                            "-species pig " .
                            "$output_dir/$idv_folder_name.filtered.mmseqs.final.fa";
    print "Start use cmd: \'$cmd_repeat_masker\'.\n";
    system $cmd_repeat_masker
        and die "Error: Command \'$cmd_repeat_masker\' failed to run normally: $?\n";


  	samtools faidx /storage3/duh/pan_genome/Anqing/$i/centrifuge/filtered.mmseqs.final.fa


    # Read the software path and set it to the default value
    my $samtools = $par->val('TOOLS', 'samtools');

    my $cmd_samtools = "$samtools " .
                       "faidx " .
                       "$output_dir/$idv_folder_name.filtered.mmseqs.final.fa";
    print "Start use cmd: \'$cmd_samtools\'.\n";
    system $cmd_samtools
        and die "Error: Command \'$cmd_samtools\' failed to run normally: $?\n";

    return 1;
}

1;

__END__
