#!/usr/bin/env perl

# Description:
# Author: zhuoy
# Date: 2021-07-21

package CDpan::Module::Mope;

use strict;
use warnings;

use Config::IniFiles;
use File::Spec::Functions qw /:ALL/;

use CDpan::Print qw / :PRINT /;

sub Mope {
    (my $par, my $idv_name) = @_;

    my $output_dir = catdir($par->val('CDPAN', 'work_dir'),'mope', $idv_name);
    mkdir $output_dir or PrintErrorMessage("Cannot create direction $output_dir: $!");

    my $output_file_prefix = catfile($output_dir, $idv_name);

    my $input_file_prefix = catfile($par->val('CDPAN', 'input_dir'), $idv_name, $idv_name);
    if ($main::modules{ "mope" }){
        unless ( -e "${input_file_prefix}.final.genome.scf.fasta"){
            PrintErrorMessage("The input file ${input_file_prefix}.final.genome.scf.fasta does not exist, whether the input direction is the output direction of Module assembly");
        }
    }

    my $min_length = $par->val('MOPE', 'min-length');
    $min_length += 0;

    PrintProcessMessage("screen sequences more than $min_length");

    # Screen more than 1k sequences
    open FASTA, '<', "${input_file_prefix}.final.genome.scf.fasta"
        or PrintErrorMessage("Couldn't open file ${input_file_prefix}.final.genome.scf.fasta: $!\n");
    open LARGE, '>', "${output_file_prefix}.final.genome.scf.large.fasta"
        or PrintErrorMessage("Couldn't create file ${output_file_prefix}.final.genome.scf.large.fasta: $!\n");

    my @fasta_new = ();
    my $fasta_length = 0;
    while (<FASTA>) {
        if (m/^[>;]/) {
            if ($fasta_length > $min_length) {
                print LARGE @fasta_new;
            }
            @fasta_new = ();
            $fasta_length = 0;
        }

        push @fasta_new, $_;
        chomp;
        $fasta_length += length;
    }
    if ($fasta_length > $min_length) {
        print LARGE @fasta_new;
    }
    close FASTA;
    close LARGE;

    # Read the software path and set it to the default value
    my $centrifuge = $par->val('TOOLS', 'centrifuge');
    my $centrifuge_kreport = $par->val('TOOLS', 'centrifuge-kreport');

    my $index = $par->val('DATA', 'index');
    my $host_taxids = $par->val('MOPE', 'host-taxids');

    my $cmd_centrifuge = "$centrifuge " .
                         "--report-file ${output_file_prefix}.centrifuge.report " .
	                     "-x $index " .
	                     "-k 1 " .
	                     "--host-taxids $host_taxids " .
	                     "-f ${output_file_prefix}.final.genome.scf.large.fasta " .
	                     "> ${output_file_prefix}.centrifuge.output " .
                         "2> ${output_file_prefix}.centrifuge.output.log";
    # print "Start use cmd: \'$cmd_centrifuge\'.\n";
    PrintProcessMessage("classifier to %%", "${output_file_prefix}.centrifuge.output");
    system $cmd_centrifuge
        and PrintErrorMessage("Command \'$cmd_centrifuge\' failed to run normally: $?\n");

    my $cmd_centrifuge_kreport = "$centrifuge_kreport " .
	                     "-x $index " .
	                     "${output_file_prefix}.centrifuge.output " .
	                     "--min-score 0 " .
	                     "--min-length 0 " .
	                     "> ${output_file_prefix}.centrifuge.krakenOut " .
                         "2> ${output_file_prefix}.centrifuge.krakenOut.log";
    # print "Start use cmd: \'$cmd_centrifuge_kreport\'.\n";
    PrintProcessMessage("make a Kraken-style report to %%", "${output_file_prefix}.centrifuge.krakenOut");
    system $cmd_centrifuge_kreport
        and PrintErrorMessage("Command \'$cmd_centrifuge_kreport\' failed to run normally: $?\n");

    if ( $par->val('CDPAN', 'output_level') == 0 ) {
        foreach (File::Slurp::read_dir($output_dir, prefix => 1)){
            if ( $_ ne "${output_file_prefix}.final.genome.scf.large.fasta" and
                $_ ne "${output_file_prefix}.centrifuge.output" and
                $_ ne "${output_file_prefix}.centrifuge.krakenOut" ){
                unlink $_ or PrintErrorMessage("Cannot delete file: $_: $!");
            }
        }
    }

    return 1;
}

1;

__END__
