#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-07-21

package CDpan::Judge;

use strict;
use warnings;
use Config::IniFiles;
use Bio::SeqIO;

sub judge {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    my $taxid = $par->val('CENTRIFUGE', 'taxid');

    open TAXID, '<', $taxid
        or PrintErrorMessage("Couldn't open file $taxid: $!\n");
    my %id_taxid = map {chomp; $_, 1} (<TAXID>);

    open INPUT, "<", "$output_dir/$idv_folder_name.centrifuge.output"
        or PrintErrorMessage("Couldn't open file $output_dir/$idv_folder_name.centrifuge.output: $!\n");
    my %id_centrifuge;
    while(<INPUT>) {
        next if m/^readID/;
        chomp;
        my @a = split m/\s+/;
        if( $id_taxid{ $a[2] } ) {
            $id_centrifuge{$a[0]} = 1;
        }
    }
    close INPUT;

    open my $INPUT, "<", "$output_dir/$idv_folder_name.final.genome.scf.large_1000.fasta"
        or PrintErrorMessage("Couldn't open file $output_dir/$idv_folder_name.final.genome.scf.large_1000.fasta: $!\n");
    open my $OUTPUT, ">", "$output_dir/$idv_folder_name.filtered.fa"
        or PrintErrorMessage("Couldn't create file $output_dir/$idv_folder_name.filtered.fa: $!\n");

    my $seq_input = Bio::SeqIO->new(-fh     => $INPUT,
                                    -format => 'Fasta');
    my $seq_output = Bio::SeqIO->new(-fh     => $OUTPUT,
                                     -format => 'Fasta');

    while ( my $seq = $seq_input->next_seq ) {
        my $seq_id = $seq->id;
        if ( $id_centrifuge{ $seq_id } ) {
            $seq_output->write_seq($seq);
        }
    }

    close $INPUT;
    close $OUTPUT;

    return 1;
}

1;

__END__
