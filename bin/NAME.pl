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
use NAME::Check;
use NAME::QualityControl;

our ($opt_h, $opt_v, $opt_d); #indicator for DEBUG mode
getopts('hvd') or abort();



# Called by Getopt::Std when supplying --help option
sub HELP_MESSAGE {
  say get_help_message();
}

# Called by Getopt::Std when supplying --version or --help option
sub VERSION_MESSAGE {
  say "Version $VERSION";
}

# Abort script due to error (e.g. invalid command line options)
sub abort {
  say get_help_message();
  exit 1;
}

# Central help message
sub get_help_message {
  return "Hi, I'm the help message :)";
}


# Abort script due to error (e.g. invalid command line options)
sub abort {
  say get_help_message();
  exit 1;
}

# Central help message
sub get_help_message {
  return "Hi, I'm the help message :)";
}
