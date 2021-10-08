#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-08-11

package CDpan::Recode;

use strict;
use warnings;

sub recode {
    # $work_dir is a folder path which contain fasqa file
    # $idv_names is a quotation of the list of the names of the individual
    (my $work_dir, my $idv_names) = @_;

    my @idv_names = @$idv_names;

    open my $OUTPUT , ">", "$work_dir/all.fasta"#TODO merge
        or die "Error: Couldn't create output file $work_dir/all.fasta: $!.\n";

    my $index = 0;
    foreach my $idv_name (@idv_names) {
        open my $INPUT, "<", "$work_dir/$idv_name.fasta"
            or die "Error: Could not open individual file $work_dir/$idv_name.fasta: $!.\n";

        while (<$INPUT>) {
            if ( m/^>/ ){
                $index += 1;
                print $OUTPUT ">cdpan$index\n";
            }else{
                print $OUTPUT $_;
            }
        }

        close $INPUT;
    }

    close $OUTPUT;

    return 1;
}

1;

__END__
