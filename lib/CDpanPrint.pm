#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-10-09

package CDpanPrint;

use strict;
use warnings;

sub PrintExitMessage {
    my $print_message = shift;

    print STDERR $print_message;
    print STDERR "\n";
    print STDERR "For more information, try \'CDpan --help\'\n";

    exit(-1);
}

1;

__END__
