#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpan::Module::Align;

use strict;
use warnings;

use Config::IniFiles;

sub Align {
    (my $par, my $idv_name) = @_;

    my $output_dir = catdir($par->val('CDPAN', 'work_dir'),'filter', $idv_name);
    mkdir $output_dir or PrintErrorMessage("Cannot create direction $output_dir: $!");

    my $output_file_prefix = catfile($output_dir, $idv_name);

    my @idv_file;
    foreach my $file (File::Slurp::read_dir( catdir($par->val('CDPAN', 'input_dir'), $idv_name ), prefix => 1)) {
        if ( $file =~ m/\S+_[12]\.fq\.gz$/ ){
            push @idv_file, $file;
        }else{
            PrintWarnMessage("Incorrectly formatted genotype data: $file, ignore it");
        }
    }



    (my $par, my $idv_folder_name, my $output_dir, my $work_dir) = @_;

    # Read the software path and set it to the default value
    my $bowtie2 = $par->val('TOOLS',  'bowtie2');
    my $thread  = $par->val('CDPAN',  'thread');

    my $cmd_align = "$bowtie2 " .
                    "-x $work_dir/index " .
                    "-U $output_dir/$idv_folder_name.singleUnmapped_R2.fq,$output_dir/$idv_folder_name.singleUnmapped_R2.fq " .
                    "-S $output_dir/$idv_folder_name.readContigAlignment.final.sam " .
                    "-p $thread " .
                    "2> $output_dir/$idv_folder_name.align.log";
    print "Start use cmd: \'$cmd_align\'.\n";
    system $cmd_align
        and PrintErrorMessage("Command \'$cmd_align\' failed to run normally: $?\n");

    return 1;
}

1;

__END__
