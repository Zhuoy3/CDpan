#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-07-21

package CDpan::Judge;

use strict;
use warnings;
use Config::IniFiles;

sub test {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    $taxid = $par->val('CENTRIFUGE', 'taxid');

    open TAXID, '<', $taxid
        or die "Error: Couldn't open file $taxid: $!\n";
    my %id_taxid;
    while(<TAXID>) {
        chomp;
        $id_taxid{$_} = "";
    }
    close TAXID;

    open INPUT, "<", "$output_dir/$idv_folder_name.centrifuge.output"
        or die "Error: Couldn't open file $output_dir/$idv_folder_name.centrifuge.output: $!\n";
    my %id_centrifuge;
    while(<INPUT>) {
        unless (m/^readID/) {
            chomp;
            my @a = split m/\s+/;
            if( exists( $id_taxid{ $a[2] } ) ) {
                $id_centrifuge{$a[0]} = "";
            }
        }
    }
    close INPUT;

    open INPUT, "<", "$output_dir/$idv_folder_name.final.genome.scf.large_1000.fasta"
        or die "Error: Couldn't open file $output_dir/$idv_folder_name.final.genome.scf.large_1000.fasta: $!\n";
    open OUTPUT, ">", "$output_dir/$idv_folder_name.filtered.fa"
        or die "Error: Couldn't create file $output_dir/$idv_folder_name.filtered.fa: $!\n";

    my $controller_output = 0;
    while(<INPUT>) {
        if( m/^[>;]/ ) {
            if( exists($id_centrifuge{ $_ }) ) {
                print OUT;
                $controller_output = 1;
            } else {
                $controller_output = 0;
            }
        }
        else {
            if($controller_output == 1) {
                print OUT;
            }
        }
    }
    close INPUT;
    close OUTPUT;

    return 1;
}

1;

__END__
