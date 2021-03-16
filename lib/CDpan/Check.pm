#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-02


package CDpan::Check;

use strict;
use warnings;
use Config::IniFiles;
use CDpan::GetPar;

sub checkpar {
    # &checkpar($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check the config file is correct
    # Because the default template is imported, the parameters will not be missing
    # Check for parameters which is redundant (may be misspelled) and missing (no default value)
    my $par = shift;

    CDpan::Check::checkparredun($par);
    CDpan::Check::checkparmissing($par);

    return 0;
}

sub checkparredun {
    # &checkpar($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check for parameters which is redundant (may be misspelled)
    my $par = shift;
    my $redundant = 0;

    my @par_sections = $par->Sections;

    foreach my $section (@par_sections) {
        if (grep { $_ eq $section } ( keys %CDpan::GetPar::par_accept )) {
            my @par_parameters = $par->Parameters($section);

            foreach my $param (@par_parameters) {
                unless (grep { $_ eq $param } @{ $CDpan::GetPar::par_accept{$section} }){
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

    return $redundant;
}

sub checkparmissing {
    # &checkpar($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check for parameters which is missing (no default value)
    my $par = shift;
    my $findmiss = 0;

    foreach my $key ( keys %CDpan::GetPar::par_require ) {
        foreach my $value ( @{ $CDpan::GetPar::par_require{$key} } ) {
            unless ( defined $par->val($key, $value)) {
                print STDERR "ERROR: Following parameters should have values:\n" unless $findmiss;
                $findmiss = 1;
                print STDERR "ERROR:     $key => $value\n";
            }
        }

    }
    die "ERROR: Please modify the parameter file and add the expected value" if $findmiss;

    return $findmiss;
}

1;

__END__
