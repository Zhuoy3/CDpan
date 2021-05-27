#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-02

package CDpan::QualityControl;

use strict;
use warnings;
use File::Spec::Functions qw /:ALL/;
use Config::IniFiles;
use File::Slurp;

sub QualityControl {
    # &qualitycontrol($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a directory path which is used to output_dir

    (my $par, my $idv_folder, my $idv_folder_name, my $output_dir) = @_;


    my @idv_file = sort ( File::Slurp::read_dir($idv_folder, prefix => 1) );
    die "ERROR: The number of genome data in the folder \'@idv_file\' is abnormal (singular)." if ( @idv_file & 1);
    foreach my $file (@idv_file) {
        die "ERROR: Incorrectly formatted genotype data: \'$file\'." unless ( $file =~ m/\S+_[12]\.fq\.gz$/ );
    }

    my $name;
    my $count = 1;
    foreach my $file (@idv_file) {
        if ( $count & 1 ) {
            die "ERROR: Two-paired data does not match: \'$file\'." unless ( ($name) = $file =~ m/\S+\/(\S+)_1\.fq\.gz$/ );
        } else {
            die "ERROR: Two-paired data does not match: \'$file\'." unless ( $file =~ m/\S+\/${name}_2\.fq\.gz$/ );
        }

       $count += 1;
    }

    my $idv_file_amount = @idv_file;

    # The default value is set when the parameter is imported, and there is no need to handle it here
    my $quality    = $par->val('QUALITYCONTROL', 'quality'); # Trim low-quality ends from reads in addition to adapter removal.
    my $length     = $par->val('QUALITYCONTROL', 'length'); # Discard reads that became shorter than length INT because of either quality or adapter trimming.
    my $error_rate = $par->val('QUALITYCONTROL', 'error_rate'); # Maximum allowed error rate (no. of errors divided by the length of the matching region)
    my $cores      = $par->val('QUALITYCONTROL', 'cores'); # Number of cores to be used for trimming

    #TODO 考虑将所有软件集中至统一模块进行调整
    # Read the software path and set it to the default value
    (my $trim_galore = $par->val('TOOLS', 'trim_galore', './trim_galore') ) =~ s/\/trim_galore$//;
    (my $cutadapt    = $par->val('TOOLS', 'cutadapt',    './cutadapt')    ) =~ s/\/cutadapt$//;
    (my $fastqc      = $par->val('TOOLS', 'fastqc',      './fastqc')      ) =~ s/\/fastqc$//;

    # Add environment variables
    $ENV{PATH} = "$trim_galore:$cutadapt:$fastqc:$ENV{PATH}:";

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
                           "--output_dir $output_dir " .
                           "-j $cores";
        #TODO 关闭trim_galore输出
        print "Start use cmd: \'$cmd_trim_galore\'\n.";
        system $cmd_trim_galore
            and die "Error: Command \'$cmd_trim_galore\' failed to run normally: $?.\n";
    }

    # merge data
    my @result_files_all = sort ( File::Slurp::read_dir($output_dir, prefix => 1) );
    my @result_files_1;
    my @result_files_2;
    foreach my $file (@result_files_all) {
        if ( $file =~ m/\Q_1_val_1.fq.gz\E$/ ) {
            push @result_files_1, $file;
        }elsif ( $file =~ m/\Q_2_val_2.fq.gz\E$/ ) {
            push @result_files_2, $file;
        }
    }
    unless ( (@result_files_1 == @result_files_2) && ( @result_files_1 == ($idv_file_amount >> 1) ) ) {
        die "ERROR: The number of result files of trim_galore is abnormal, trim_galore may not be running normally.";
    }

    @result_files_1 = sort @result_files_1;
    @result_files_2 = sort @result_files_2;

    my $cmd_merge_data_1 = "cat @result_files_1 > $output_dir/${idv_folder_name}_clean_1.fq.gz";
    print "Start use cmd: \'$cmd_merge_data_1\'.\n";
    system $cmd_merge_data_1
        and die "Error: Command \'$cmd_merge_data_1\' failed to run normally: $?.\n";
    unlink @result_files_1;

    my $cmd_merge_data_2 = "cat @result_files_2 > $output_dir/${idv_folder_name}_clean_2.fq.gz";
    print "Start use cmd: \'$cmd_merge_data_2\'.\n";
    system $cmd_merge_data_2
        and die "Error: Command \'$cmd_merge_data_2\' failed to run normally: $?.\n";
    unlink @result_files_2;

    return 1;
}


1;

__END__
