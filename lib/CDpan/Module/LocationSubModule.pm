#!/usr/bin/env perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpan::Module::LocationSubModule;

use strict;
use warnings;

use Config::IniFiles;
use File::Spec::Functions qw /:ALL/;

use CDpan::Print qw / :PRINT /;

sub RepeatMasker {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    (my $par) = @_;

    my $output_dir = catdir($par->val('CDPAN', 'work_dir'), 'repeat_masker');

    my $input_file = catfile($par->val('CDPAN', 'input_dir'), 'merge', 'merge.fasta');
    if ($main::modules{ "location" }){
        unless ( -e $input_file){
            PrintErrorMessage("The input file $input_file does not exist, whether the input direction is the output direction of Module merge");
        }
    }

    copy $input_file, "$output_dir/merge.fasta"
        or PrintErrorMessage("Cannot copy file $input_file to $output_dir/merge.fasta: $!");

    # Read the software path and set it to the default value
    my $repeat_masker = $par->val('TOOLS', 'RepeatMasker');

    my $thread  = $par->val('CDPAN',    'thread');
    my $species = $par->val('LOCATION', 'species');

    my $cmd_repeat_masker = "$repeat_masker " .
                            "-parallel $thread " .
                            "-nolow " .
                            "-species $species " .
                            "$output_dir/merge.fasta " .
                            "> $output_dir/repeat_masker.log";
    print "Start use cmd: \'$cmd_repeat_masker\'.\n";
    system $cmd_repeat_masker
        and PrintErrorMessage("Command \'$cmd_repeat_masker\' failed to run normally: $?\n");

    # Read the software path and set it to the default value
    my $samtools = $par->val('TOOLS', 'samtools');

    my $cmd_samtools = "$samtools " .
                       "faidx " .
                       "$output_dir/merge.fasta";
    print "Start use cmd: \'$cmd_samtools\'.\n";
    system $cmd_samtools
        and PrintErrorMessage("Command \'$cmd_samtools\' failed to run normally: $?\n");

    # Read the software path and set it to the default value
    my $bowtie2_build = $par->val('TOOLS', 'bowtie2-build');
    my $cmd_bowtie2_build = "$bowtie2_build " .
                            "$output_dir/merge.fasta.masked " .
                            "$output_dir/index " .
                            "> $output_dir/bowtie2_build1.log " .
                            "2> $output_dir/bowtie2_build2.log";
    print "Start use cmd: \'$cmd_bowtie2_build\'.\n";
    system $cmd_bowtie2_build
        and PrintErrorMessage("Command \'$cmd_bowtie2_build\' failed to run normally: $?\n");

    return 1;
}

1;

__END__
