#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpan;

use strict;
use warnings;

use CDpan::Print qw \ :NONE \ ;
use CDpan::PreProcess;
use CDpan::Master;


require Exporter;
our @ISA = qw \ Exporter \;
our @EXPORT = qw \ PrintExitMessage PrintWarnMessage PrintErrorMessage PrintStartMessage PrintEndMessage
                   PreProcess
                   Filter Align Extract Assembly Mope Vot Soot Merge Location RunAll RunDisplace\;
our %EXPORT_TAGS = (
    ALL    => [ @EXPORT ],
    PRINT  => [ qw \ PrintExitMessage PrintWarnMessage PrintErrorMessage PrintStartMessage PrintEndMessage\ ],
    MODULE => [ qw \ Filter Align Extract Assembly Mope Vot Soot Merge Location RunAll RunDisplace\ ],
);

# module CDpan::Print
sub PrintExitMessage  { return &CDpan::Print::PrintExitMessage  };
sub PrintWarnMessage  { return &CDpan::Print::PrintWarnMessage  };
sub PrintErrorMessage { return &CDpan::Print::PrintErrorMessage };
sub PrintStartMessage { return &CDpan::Print::PrintStartMessage };
sub PrintEndMessage   { return &CDpan::Print::PrintEndMessage   };

# module CDpan::Check
sub PreProcess { return &CDpan::PreProcess::PreProcess };

# module CDpan::Master
sub Filter      { return &CDpan::Master::Filter    };
sub Align       { return &CDpan::Master::Align     };
sub Extract     { return &CDpan::Master::Extract   };
sub Assembly    { return &CDpan::Master::Assembly  };
sub Mope        { return &CDpan::Master::Mope      };
sub Vot         { return &CDpan::Master::Vot       };
sub Soot        { return &CDpan::Master::Soot      };
sub Merge       { return &CDpan::Master::Merge     };
sub Location    { return &CDpan::Master::Location  };
sub RunAll      { return &CDpan::Master::RunAll    };
sub RunDisplace { return &CDpan::Master::RunDisplace };

1;

__END__
