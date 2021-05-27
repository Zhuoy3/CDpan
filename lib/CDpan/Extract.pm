#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-17

package CDpan::Extract;

use strict;
use warnings;
use Config::IniFiles;

sub extract {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a directory path which is used to output

    (my $par, my $idv_folder_name, my $output_dir) = @_;

    #TODO 考虑将所有软件集中至统一模块进行调整
    # Read the software path and set it to the default value
    (my $gatk = $par->val('TOOLS', 'gatk', './gatk') ) =~ s/\/gatk$//;
    $ENV{PATH} = "$gatk:$ENV{PATH}:";

    my $ref = $par->val('DATA', 'ref');

    unless ($main::ref_dict) {
        if ( -e "$ref.dict" ) {
            rename "$ref.dict" => "$ref.dict.old" or die "ERROR: Cannot change the name of file '$ref.dict': $!\n";
            warn "WARNING: File '$ref.dict' exists and has been renamed to '$ref.dict.old'";
        }
        my $cmd_gatk_dict = "gatk --java-options \"-Xmx32G\" CreateSequenceDictionary " .
                            "-R $ref " .
                            "-O $ref.dict";
        print "Start use cmd: \'$cmd_gatk_dict\'.\n";
        system $cmd_gatk_dict
            and die "Error: Command \'$cmd_gatk_dict\' failed to run normally: $?.\n";
        $main::ref_dict = 1;
    }

    #reorder sam
    my $cmd_reorder = "gatk --java-options \"-Xmx32G\" ReorderSam " .
                      "-I $output_dir/$idv_folder_name.sam " .
                      "-O $output_dir/$idv_folder_name.reorder.sam " .
                      "-R $ref " .
                      "-SD $ref.dict";
    print "Start use cmd: \'$cmd_reorder\'.\n";
    system $cmd_reorder
        and die "Error: Command \'$cmd_reorder\' failed to run normally: $?.\n";

    #TODO 考虑将所有软件集中至统一模块进行调整
    # Read the software path and set it to the default value
    (my $samtools = $par->val('TOOLS', 'samtools', './samtools') ) =~ s/\/samtools$//;
    $ENV{PATH} = "$samtools:$ENV{PATH}:";

    #sam to bam
    my $cmd_sam2bam = "samtools view -b  " .
                      "-S $output_dir/$idv_folder_name.reorder.sam " .
                      "-o $output_dir/$idv_folder_name.reorder.bam";
    print "Start use cmd: \'$cmd_sam2bam\'.\n";
    system $cmd_sam2bam
        and die "Error: Command \'$cmd_sam2bam\' failed to run normally: $?.\n";
    unlink "$output_dir/$idv_folder_name.reorder.sam";

    #sort bam
    my $cmd_sort = "gatk --java-options \"-Xmx32G\" SortSam " .
                   "-I $output_dir/$idv_folder_name.reorder.bam " .
                   "-O $output_dir/$idv_folder_name.bam " .
                   "--SORT_ORDER coordinate";
    print "Start use cmd: \'$cmd_sort\'.\n";
    system $cmd_sort
        and die "Error: Command \'$cmd_sort\' failed to run normally: $?.\n";
    unlink "$output_dir/$idv_folder_name.reorder.bam";

    return 1;
}


1;

__END__
