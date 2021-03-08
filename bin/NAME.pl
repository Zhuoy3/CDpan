#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-03

use strict;
use warnings;
use feature qw(say);

use FindBin;
use lib "$FindBin::Bin/../lib";

use Carp; #cluck(warn) confess(die)
use Getopt::Std;
#use NAME::Check;
#use NAME::QualityControl;

our $VERSION = '0.0.1';

our ($opt_h, $opt_v, $opt_d); #indicator for DEBUG mode
getopts('hvd') ;#or abort();



# Called by Getopt::Std when supplying --help option
sub HELP_MESSAGE {
    system "pod2text $FindBin::Bin/../doc/help.pod";
    exit 0;
}

__END__

# Called by Getopt::Std when supplying --version
sub VERSION_MESSAGE {
    say "Version $VERSION"; # !! maybe need other thing
    exit 0;
}

__END__

# Abort script due to error
sub abort {
  say get_help_message();
  exit 1;
}

# Central help message
sub get_message {
    my $file_path = shift;
    my $file = `pod2text $file_path`;
    return "Hi, I'm the help message :)";
}
