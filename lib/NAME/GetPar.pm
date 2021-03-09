#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-07

package NAME::GetPar;

use strict;
use warnings;
use Config::IniFiles;

sub getpar {
    # &getpar($opt)
    my $file_par_path = shift;

    # Since "FindBin" uses a "BEGIN" block, it'll be executed only once, and only the first caller will get it right.
    # $FindBin::Bin is path to bin directory from where 'NAME.pl' was invoked
    my $par_default;
    $par_default = new Config::IniFiles(-file => "$FindBin::Bin/../config/default.ini");

    my $par;
    $par = new Config::IniFiles(-file                => $file_par_path,
                                -import              => $par_default,
                                -allowedcommentchars => '#')
        or die("ERROR: Could not import parameter from $file_par_path: @Config::IniFiles::errors");

    return $par;
}

sub printpar {
    # use in debug mode

}

1;

__END__
