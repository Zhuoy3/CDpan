#!/usr/bin/env perl

# Description:
# Author: zhuoy
# Date: 2021-10-09

package CDpan::Print;

use strict;
use warnings;

use File::Spec::Functions qw /:ALL/;
use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

require Exporter;
our @ISA = qw \ Exporter \;
our @EXPORT = qw \ PrintExitMessage PrintWarnMessage PrintErrorMessage PrintStartMessage PrintEndMessage PrintProcessMessage \;
our %EXPORT_TAGS = (
     PRINT => [ qw \ PrintExitMessage PrintWarnMessage PrintErrorMessage PrintStartMessage PrintEndMessage PrintProcessMessage\ ],
     NONE => []);

sub PrintExitMessage {
    my $print_message = shift;

    print STDERR "$print_message\n";
    print STDERR "\n";
    print STDERR "For more information, try \'CDpan --help\'\n";

    exit(-1);
}

sub PrintWarnMessage {
    my $print_message = shift;

    print STDERR BOLD MAGENTA "Warning: ";
    print STDERR "$print_message\n";

    return 1;
}

sub PrintErrorMessage {
    my $print_message = shift;
    my $is_die = shift // 1;

    print STDERR "\n";
    print STDERR BOLD RED "Error: ";
    print STDERR "$print_message\n";
    print STDERR "\n";

    exit(255) if $is_die;

    return 1;
}

sub PrintStartMessage {
    my $print_message = shift;

    print STDERR "--------------------------------------------------------------------------------\n";
    print STDERR "\n";
    print STDERR "$print_message\n";
    print STDERR "\n";

    return 1;
}

sub PrintEndMessage {
    my $print_message = shift;

    print STDERR "\n";
    print STDERR "$print_message\n";
    print STDERR "\n";

    return 1;
}

sub PrintProcessMessage {
    my @print_message = split /%%/, shift;
    my @print_files = @_;

    print STDERR "    ";

    if ( @print_message == 1 and @print_files == 0) {
        print STDERR "@print_message";
    }
    else{
        foreach (@print_message) {
            print STDERR $_;

            my $file = shift @print_files;

            if ( ref $file ){
                foreach my $file_sub (@{$file}){
                    ( my $file_name, my $idv_dir_name ) = reverse splitdir($file_sub);
                    print STDERR "$idv_dir_name/$file_name ";
                }
            }
            else {
                ( my $file_name, my $idv_dir_name ) = reverse splitdir($file);
                print STDERR "$idv_dir_name/$file_name";
            }
        }
    }

    print STDERR "\n";

    return 1;
}



1;

__END__
