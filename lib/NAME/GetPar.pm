#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-07

package NAME::GetPar;

use strict;
use warnings;
use Data::Dumper;
use Config::IniFiles;

sub getpar {
    # &getpar($opt)
    # Import parameter form file which be given by command line
    my $file_par_path = shift;

    # Since "FindBin" uses a "BEGIN" block, it'll be executed only once, and only the first caller will get it right
    # $FindBin::Bin is path to bin directory from where 'NAME.pl' was invoked
    my $par_default;
    $par_default = new Config::IniFiles(-file => "$FindBin::Bin/../config/default.ini");

    my $par;
    $par = new Config::IniFiles(-file                => $file_par_path,
                                -import              => $par_default,
                                -allowedcommentchars => '#')
        or die "ERROR: Could not import parameter from $file_par_path: @Config::IniFiles::errors";

    NAME::GetPar::checkpar($par, $par_default);

    $par->WriteConfig("$file_par_path.imp") if $main::opt_d;

    return $par;
}

sub checkpar {
    # &checkpar($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check the config file is correct
    # Because the default template is imported, the parameters will not be missing
    # Check for parameters which is redundant (may be misspelled) and no-value (no default value)
    my $par = shift;
    my $redundant = 0;

    my @par_sections = $par->Sections;

    foreach my $section (@par_sections) {
        if (grep { $_ eq $section } ( keys %NAME::GetPar::needpar )) {
            my @par_parameters = $par->Parameters($section);

            foreach my $param (@par_parameters) {
                unless (grep { $_ eq $param } @{ $NAME::GetPar::needpar{$section} }){
                    warn "WARNING: Unknown parameter found: '[$section] => $param', will be removed.\n";
                    $par->delval($section, $param);
                    die "ERROR: More unknown parameter, please check the parameter file.\n" if (++$redundant > 3);
                }
            }
        }else{
            warn "WARNING: Unknown parameter found: '[$section] => ...', will be removed.\n";
            $par->DeleteSection($section);
            die "ERROR: More unknown parameter, please check the parameter file.\n" if (++$redundant > 3);
        }
    }

    my $findmiss = 0;
    foreach my $key ( keys %NAME::GetPar::needvalues ) {
        foreach my $value ( @{ $NAME::GetPar::needvalues{$key} } ) {
            if ( $par->val($key, $value) eq 'null') {
                print STDERR "ERROR: Following parameters should have values:\n" unless $findmiss;
                $findmiss = 1;
                print STDERR "ERROR:     $key => $value\n";
            }
        }

    }
    die "ERROR: Please modify the parameter file and add the expected value" if $findmiss;

    return 0;
}

our %needpar = (
    'DATA' => ['a'],
    'QC' => ['t'],
    'TRIM' => ['b'],
);

our %needvalues = (
    'QC' => ['t'],
);

1;

__END__
