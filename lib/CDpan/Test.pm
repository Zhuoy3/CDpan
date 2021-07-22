#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
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

    # Read the software path and set it to the default value
    my $centrifuge = $par->val('TOOLS', 'centrifuge');
    my $centrifuge_kreport = $par->val('TOOLS', 'centrifuge-kreport');

    my $cmd_centrifuge = "$centrifuge " .
                         "--report-file $output_dir/$idv_folder_name.centrifuge.report " .
	                     "-x /diskd/duh/domestic/pan_genome/centrifuge/index/nt " .
	                     "-k 1 " .
	                     "--host-taxids 9823 " .
	                     "-f $output_dir/$idv_folder_name.final.genome.scf.fasta " .
	                     "> $output_dir/$idv_folder_name.centrifuge.output";
    #TODO where is -x???
    print "Start use cmd: \'$cmd_centrifuge\'.\n";
    system $cmd_centrifuge
        and die "Error: Command \'$cmd_centrifuge\' failed to run normally: $?.\n";

    my $centrifuge_kreport = "$centrifuge_kreport " .
	                     "-x /diskd/duh/domestic/pan_genome/centrifuge/index/nt " .
	                     "$output_dir/$idv_folder_name.centrifuge.output " .
	                     "--min-score 0 " .
	                     "--min-length 0 " .
	                     "> $output_dir/$idv_folder_name.centrifuge.krakenOut";
    #TODO where is -x???
    print "Start use cmd: \'$cmd_centrifuge\'.\n";
    system $cmd_centrifuge
        and die "Error: Command \'$cmd_centrifuge\' failed to run normally: $?.\n";

    return 1;
}

1;

__END__
