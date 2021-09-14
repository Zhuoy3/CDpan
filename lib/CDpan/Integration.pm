#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-08-21

package CDpan::Integration;

use strict;
use warnings;
use Config::IniFiles;
use File::Copy qw / copy /;

sub integration {
    # &extract($opt, $idv_folder_name, $output_dir, $work_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir, my $work_dir) = @_;

    mkdir "$output_dir/$link_new"
        or die "Error: Cannot create process folder '$output_dir/$link_new': $!\n";

    # switch to destination file
    copy "$output_dir/$idv_folder_name.mateLinks.txt", "$output_dir/$link_new/mateLinks.txt"
        or die "Error: Cannot copy file '$output_dir/$idv_folder_name.mateLinks.txt': $!\n";
    copy "$work_dir/all.fasta.fai", "$output_dir/$link_new/all.fasta.fai"
        or die "Error: Cannot copy file '$work_dir/all.fasta.fai': $!\n";

    chdir "$output_dir/$link_new"
        or die "Error: Cannot chdir to '$output_dir/$link_new': $!\n";

	mkdir "./1";
	mkdir "./2";
	mkdir "./3";
	mkdir "./4";

    # archive the contig name and length
    system "awk '{print \$1,\$2}' ./all.fasta.fai > ./contig.name"
        and die "Error: Command failed to run normally: $?\n";

    # 1do
    open my $CONTIG, '<', "./mateLinks.txt"
        or die "Error: Cannot open file '$output_dir/$link_new/contig.name': $!\n";

    my %contig;
    my $line = 0;
    while(<$CONTIG>) {
        chomp;
        my @contig_lines = split /\s+/, $_;
        $contig{$contig_lines[1]}{$line} = $_;
        $line++;
    }
    close $CONTIG;


    open my $INPUT, "<", "./contig.name"
        or die "Error: Cannot open file '$output_dir/$link_new/mateLinks.txt': $!\n";
    open my $OUTPUT1, ">", "./1/1.name"
        or die "Error: Cannot create file '$output_dir/$link_new/1/1.name': $!\n";
    open my $OUTPUT2, ">", "./1/2to4.name"
        or die "Error: Cannot create file '$output_dir/$link_new/1/2to4.name': $!\n";

    while(<$INPUT>) {
        chomp;
        my @lines_links = split /\s+/, $_;
        my $have_contigs = 0;

        open my $LINK, ">", "./4/$lines_links[0].link"
            or die "Error: Cannot create file '$output_dir/$link_new/4/$lines_links[0].link': $!\n";

        if( exists $contig{$lines_links[0]} ) {
            foreach my $key (keys %{$contig{$lines_links[0]}}) {
                print $LINK "$contig{$lines_links[0]}{$key}\n";
            }
            $have_contigs = 1;
        }
        close $LINK;

        if($have_contigs) {
            print $OUTPUT2 "$_\n";
        }else{
            print $OUTPUT1 "$_\n";
        }
    }

    close $INPUT;

    return 1;
}

1;

__END__
