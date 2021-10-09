#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-08-25

package CDpan::DeRepeat;

use strict;
use warnings;
use Config::IniFiles;
use File::Copy qw / copy /;

sub de_repeat {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    if ( -s "$output_dir/$idv_folder_name.filtered.mmseqs.delta.coords" ) {
        my $cmd_repeat = "tail -n +6 $output_dir/$idv_folder_name.filtered.mmseqs.delta.coords " .
                        "| awk \'{print \$19}\' - " .
                        "| sort -u - " .
                        "> $output_dir/$idv_folder_name.repeat.names";
        print "Start use cmd: \'$cmd_repeat\'.\n";
        system $cmd_repeat
            and die "Error: Command \'$cmd_repeat\' failed to run normally: $?\n";

        open my $REPEAT, "<", "$output_dir/$idv_folder_name.repeat.names"
            or die "Error: Couldn't open file $output_dir/$idv_folder_name.repeat.names: $!\n";
        my %id_repeat = map {chomp; $_, 1} (<$REPEAT>);

        open my $INPUT, "<", "$output_dir/$idv_folder_name.filtered.mmseqs_rep_seq.fasta"
            or die "Error: Couldn't open file $output_dir/$idv_folder_name.filtered.mmseqs_rep_seq.fasta: $!\n";
        open my $OUTPUT, ">", "$output_dir/$idv_folder_name.filtered.mmseqs.final.fa"
            or die "Error: Couldn't create file $output_dir/$idv_folder_name.filtered.mmseqs.final.fa: $!\n";

        my $seq_input = Bio::SeqIO->new(-fh     => $INPUT,
                                        -format => 'Fasta');

        while (my $seq = $seq_input->next_seq) {
            my $seq_id = $seq->id;
            unless ( $id_repeat{ $seq_id } ) {
                print $OUTPUT ">$seq_id \n" . $seq->seq . "\n";
            }
        }
    }
    else{
        copy("$output_dir/$idv_folder_name.filtered.mmseqs_rep_seq.fasta", "$output_dir/$idv_folder_name.filtered.mmseqs.final.fa")
             or die "Error: Copy file \'$output_dir/$idv_folder_name.filtered.mmseqs_rep_seq.fasta\' to \'$output_dir/$idv_folder_name.filtered.mmseqs.final.fa\' failed: $!\n";
    }

    return 1;
}

1;

__END__
