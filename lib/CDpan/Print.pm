#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-10-09

package CDpanPrint;

use strict;
use warnings;

use Term::ANSIColor qw(:constants);
$Term::ANSIColor::AUTORESET = 1;

require Exporter;
our @ISA = qw \ Exporter \;
our @EXPORT = qw \ PrintExitMessage PrintWarnMessage PrintErrorMessage \;
our %EXPORT_TAGS = ( ALL => [ @EXPORT ] );

sub PrintExitMessage {
    my $print_message = shift;

    print STDERR $print_message;
    print STDERR "\n";
    print STDERR "For more information, try \'CDpan --help\'\n";

    exit(-1);
}

sub PrintWarnMessage {
    my $print_message = shift;

    print STDERR BOLD MAGENTA "Warning: ";
    print STDERR $print_message;

    return 1;
}

sub PrintErrorMessage {
    my $print_message = shift;
    my $is_die = shift // 1;

    print STDERR "\n";
    print STDERR BOLD RED "Error: ";
    print STDERR $print_message;
    print STDERR "\n";

    exit(255) if $is_die;

    return 1;
}

1;

__END__
