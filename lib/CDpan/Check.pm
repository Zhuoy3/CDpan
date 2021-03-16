#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-02


package CDpan::Check;

use strict;
use warnings;
use File::Spec::Functions;
use Config::IniFiles;
use CDpan::GetPar;

sub checkpar {
    # &checkpar($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check the config file is correct
    # Because the default template is imported, the parameters will not be missing
    # Check for parameters which is redundant (may be misspelled) and missing (no default value)
    my $par = shift;

    # !!!No parameter check for easy debugging
    CDpan::Check::checkpar_redun($par) unless $main::opt_d;
    CDpan::Check::checkpar_missing($par) unless $main::opt_d;
    CDpan::Check::checkpar_tool($par);

    return 0;
}

sub checkpar_redun {
    # &checkpar_redun($opt)
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

sub checkpar_missing {
    # &checkpar_missing($opt)
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

sub checkpar_tool {
    # &checkpar_tool($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check whether the tools specified in par is available, and search tools in PATH if not specified
    my $par = shift;
    my @tools_needed = qw \ trim_galores
                            cutadapt
                            fastqc \;
    my @tools_missing;

    FINDTOOL: foreach my $tools (@tools_needed) {
        if ( defined $par->val('TOOLS', $tools) ) {
            ( -e -x $par->val('TOOLS', $tools) ) ? next : $par->delval('TOOLS', $tools);
        }

        foreach my $path ( path() )  {
            next FINDTOOL if (-e -x catfile($path, $tools));
        }

        push @tools_missing, $tools;
    }

    die "ERROR: can't locate '@tools_missing', please add these to PATH or use 'TOOLS' section of parameters file.\n If the above is indeed executed, please confirm whether you have execution authority.\n" if (@tools_missing);

    return @tools_missing;
}

1;

__END__
