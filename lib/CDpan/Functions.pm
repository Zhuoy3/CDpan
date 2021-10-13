#!/usr/bin/env perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpan::Functions;

use strict;
use warnings;

use CDpan::Print;

require Exporter;
our @ISA = qw \ Exporter \;
our @EXPORT = qw \ PrintExitMessage PrintWarnMessage PrintErrorMessage
                   PreProcess \;
our %EXPORT_TAGS = ( ALL => [ @EXPORT ],
                     PRINT => [ qw \ PrintExitMessage PrintWarnMessage PrintErrorMessage\ ] );

sub PrintExitMessage  { return &CDpan::Print::PrintExitMessage };
sub PrintWarnMessage  { return &CDpan::Print::PrintWarnMessage };
sub PrintErrorMessage { return &CDpan::Print::PrintErrorMessage };

1;

__END__
