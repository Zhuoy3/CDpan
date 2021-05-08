#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-05-08
# called by: CDpan.pl


package CDpan::GetSample;

use strict;
use warnings;
use File::Spec::Functions;
use Config::IniFiles;
use File::Slurp;

sub GetSampleFolder {
    # &GetSampleFolder($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # $opt -> [DATA]:input should be used which is a directory path containing genomic data
    # directory contains multiple subfolders, each folder contains an individual information
    my $opt = shift;
    my $directory =  $opt->val('DATA', 'input');

    die "ERROR: Data folder \'$directory\' does not exist." unless ( -e $directory );

    # Get individual subfolders that need to be processed
    my @idv_folder = sort ( File::Slurp::read_dir($directory, prefix => 1) );

    return @idv_folder;
}

sub GetSampleFile {
    # &GetSampleFile($directory)
    # $directory is a directory path containing genomic data of a individual
    my $directory = shift;

    my @gene_folder = sort ( File::Slurp::read_dir($directory, prefix => 1) );

    die "ERROR: The number of genome data in the folder \'$directory\' is abnormal (singular)." if ( @gene_folder & 1);
    foreach my $file (@gene_folder) {
        die "ERROR: Incorrectly formatted genotype data: \'$file\'." unless ( $file =~ m/\S+_[12]\.fq\.gz$/ );
    }

    my $name;
    my $count = 1;
    foreach my $file (@gene_folder) {
        if ( $count & 1 ) {
            die "ERROR: Two-paired data does not match: \'$file\'." unless ( ($name) = $file =~ m/\S+\/(\S+)_1\.fq\.gz$/ );
        } else {
            die "ERROR: Two-paired data does not match: \'$file\'." unless ( $file =~ m/\S+\/${name}_2\.fq\.gz$/ );
        }

       $count += 1;
    }

    return @gene_folder;
}

1;

__END__
