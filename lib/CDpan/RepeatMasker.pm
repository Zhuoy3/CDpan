#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpan::RepeatMasker;

use strict;
use warnings;
use Config::IniFiles;

sub repeat_masker {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    (my $par, my $work_dir) = @_;

    # Read the software path and set it to the default value
    # my $repeat_masker = $par->val('TOOLS', 'RepeatMasker');

    # my $cmd_repeat_masker = "$repeat_masker " .
    #                         "-parallel 12 " .
    #                         "-nolow " .
    #                         "-species pig " . #TODO
    #                         "$work_dir/all.fasta " .
    #                         "> $work_dir/repeat_masker.log";
    # print "Start use cmd: \'$cmd_repeat_masker\'.\n";
    # system $cmd_repeat_masker
    #     and die "Error: Command \'$cmd_repeat_masker\' failed to run normally: $?\n";

    # Read the software path and set it to the default value
    my $samtools = $par->val('TOOLS', 'samtools');

    my $cmd_samtools = "$samtools " .
                       "faidx " .
                       "$work_dir/all.fasta";
    print "Start use cmd: \'$cmd_samtools\'.\n";
    system $cmd_samtools
        and die "Error: Command \'$cmd_samtools\' failed to run normally: $?\n";

    # Read the software path and set it to the default value
    my $bowtie2_build = $par->val('TOOLS', 'bowtie2-build');
    my $cmd_bowtie2_build = "$bowtie2_build " .
                            "$work_dir/all.fasta.masked " .
                            "$work_dir/index " .
                            "> $work_dir/bowtie2_build1.log " .
                            "2> $work_dir/bowtie2_build2.log";
    print "Start use cmd: \'$cmd_bowtie2_build\'.\n";
    system $cmd_bowtie2_build
        and die "Error: Command \'$cmd_bowtie2_build\' failed to run normally: $?\n";

    return 1;
}

1;

__END__
