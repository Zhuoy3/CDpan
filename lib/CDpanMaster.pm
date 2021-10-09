#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpanMaster;

use strict;
use warnings;

# use CDpan::Parameter;
# use CDpan::QualityControl;
# use CDpan::Comparison;
# use CDpan::Extract;
# use CDpan::Assembly;
# use CDpan::Test;
# use CDpan::Judge;
# use CDpan::MMSeqs;
# use CDpan::Nucmer;
# use CDpan::DeRepeat;
# use CDpan::Recode;
# use CDpan::RepeatMasker;
# use CDpan::Align;
# use CDpan::Change;
# use CDpan::Integration;

require Exporter;
our @ISA = qw \ Exporter \;

our @EXPORT = qw \ PreCheck \;

our %EXPORT_TAGS = ( ALL => [ @EXPORT ] );

sub PreCheck {


    return 1;
}

1;

__END__
