#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-07-21

package CDpan::Test;

use strict;
use warnings;
use Config::IniFiles;

sub test {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    # Screen more than 1k sequences
    open FASTA, '<', "$output_dir/$idv_folder_name.final.genome.scf.fasta"
        or PrintErrorMessage("Couldn't open file $output_dir/$idv_folder_name.final.genome.scf.fasta: $!\n");
    open LARGE, '>', "$output_dir/$idv_folder_name.final.genome.scf.large_1000.fasta"
        or PrintErrorMessage("Couldn't create file $output_dir/$idv_folder_name.final.genome.scf.large_1000.fasta: $!\n");

    my @fasta_new = ();
    my $fasta_length = 0;
    while (<FASTA>) {
        if (m/^[>;]/) {
            if ($fasta_length > 1000) {
                print LARGE @fasta_new;
            }
            @fasta_new = ();
            $fasta_length = 0;
        }

        push @fasta_new, $_;
        chomp;
        $fasta_length += length;
    }
    if ($fasta_length > 1000) {
        print LARGE @fasta_new;
    }
    close FASTA;
    close LARGE;

    # Read the software path and set it to the default value
    my $centrifuge = $par->val('TOOLS', 'centrifuge');
    my $centrifuge_kreport = $par->val('TOOLS', 'centrifuge-kreport');

    my $index = $par->val('CENTRIFUGE', 'index');

    my $cmd_centrifuge = "$centrifuge " .
                         "--report-file $output_dir/$idv_folder_name.centrifuge.report " .
	                     "-x $index " .
	                     "-k 1 " .
	                     "--host-taxids 9823 " .
	                     "-f $output_dir/$idv_folder_name.final.genome.scf.large_1000.fasta " .
	                     "> $output_dir/$idv_folder_name.centrifuge.output " .
                         "2> $output_dir/$idv_folder_name.centrifuge.output.log";
    print "Start use cmd: \'$cmd_centrifuge\'.\n";
    system $cmd_centrifuge
        and PrintErrorMessage("Command \'$cmd_centrifuge\' failed to run normally: $?\n");

    my $cmd_centrifuge_kreport = "$centrifuge_kreport " .
	                     "-x $index " .
	                     "$output_dir/$idv_folder_name.centrifuge.output " .
	                     "--min-score 0 " .
	                     "--min-length 0 " .
	                     "> $output_dir/$idv_folder_name.centrifuge.krakenOut " .
                         "2> $output_dir/$idv_folder_name.centrifuge.krakenOut.log";
    print "Start use cmd: \'$cmd_centrifuge_kreport\'.\n";
    system $cmd_centrifuge_kreport
        and PrintErrorMessage("Command \'$cmd_centrifuge_kreport\' failed to run normally: $?\n");

    return 1;
}

1;

__END__
