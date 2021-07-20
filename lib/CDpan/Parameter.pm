#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-07
# called by: CDpan.pl

package CDpan::Parameter;

use strict;
use warnings;
use Config::IniFiles;

sub InputPar {
    # &InputPar($opt)
    # Import parameter form file which be given by command line
    my $file_par_path = shift;

    # Since "FindBin" uses a "BEGIN" block, it'll be executed only once, and only the first caller will get it right
    # $FindBin::Bin is path to bin directory from where 'CDpan.pl' was invoked
    my  $par_default = new Config::IniFiles(-file => "$FindBin::Bin/../config/par_default.ini");
    my $par = new Config::IniFiles(-file                => $file_par_path,
                                   -import              => $par_default,
                                   -allowedcommentchars => '#')
        or die "ERROR: Could not import parameter from \'$file_par_path\': @Config::IniFiles::errors.\n";



    _CheckRedundant($par);
    _CheckMissing($par);

    #CDpan::Check::checkpar_tool($par);

    # unless ( -e $par->val('COMPARISON', 'ref') ) {
    #     my $reference = $par->val('COMPARISON', 'ref');
    #     die "ERROR: Reference genome file \'$reference\' does not exist.\n";
    # }



    #TODO 需要一个检查转换为全局路径的函数

    $par->WriteConfig("$file_par_path.import") if $main::debug;

    return $par;
}

sub _GetStdPar {
    # &_GetStdPar($file)
    # Import parameter form file in the folder '/CDpan/config/'
    my $file_path = shift;
    my %output;

    open PARFILE, "<", $file_path or die "ERROR: Cannot open file \'$file_path\': $!";
    while (<PARFILE>) {
        next unless s/^(.*?):\s*//; # Capturing the first parameter of the line: Section
        $output{$1} =[ split ]; # A hash that points to the list
    }
    close PARFILE;

    return %output;
}

sub _CheckRedundant {
    # &_CheckRedundant($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check for parameters which is redundant (may be misspelled)
    my $par = shift;
    my $redundant = 0;

    my %par_accept  = _GetStdPar("$FindBin::Bin/../config/par_accept.txt");
    my @par_sections = $par->Sections;

    foreach my $section (@par_sections) {
        if (grep { $_ eq $section } ( keys %par_accept )) {
            my @par_parameters = $par->Parameters($section);

            foreach my $param (@par_parameters) {
                unless (grep { $_ eq $param } @{ $par_accept{$section} }){
                    warn "WARNING: Unknown parameter found: '[$section] => $param', will be removed.\n";
                    $par->delval($section, $param);
                    $redundant = 1;
                }
            }
        }else{
            warn "WARNING: Unknown parameter found: '[$section] => ...', will be removed.\n";
            $par->DeleteSection($section);
            $redundant = 1;
        }
    }

    die "ERROR: More unknown parameter, please check the parameter file.\n" if $redundant;

    return $redundant;
}

sub _CheckMissing {
    # &_CheckMissing($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check for parameters which is missing (no default value)
    my $par = shift;
    my $findmiss = 0;

    my %par_require = _GetStdPar("$FindBin::Bin/../config/par_require.txt");
    foreach my $key ( keys %par_require ) {
        foreach my $value ( @{ $par_require{$key} } ) {
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

sub _CheckTools {
    # &_CheckTools($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check whether the tools specified in par is available, and search tools in PATH if not specified
    my $par = shift;

    my @tools_needed = qw \ trim_galore cutadapt fastqc bwa gatk samtools masurca \;
    my @tools_missing;

    FINDTOOL: foreach my $tools (@tools_needed) {
        if ( defined $par->val('TOOLS', $tools) ) {
            ( -e -x $par->val('TOOLS', $tools) ) ? next : $par->delval('TOOLS', $tools);
        }

        foreach my $path ( path() )  {
            if (-e -x catfile($path, $tools)) {

                next FINDTOOL;
            };
        }

        push @tools_missing, $tools;
    }

    die "ERROR: can't locate '@tools_missing', please add these to PATH or use 'TOOLS' section of parameters file.\n If the above is indeed executed, please confirm whether you have execution authority.\n" if (@tools_missing);

    return @tools_missing;
}

sub checkpar_file {
    # &checkpar_tool($opt, $file)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check if the required file exists
    my $par = shift;
    my $file = shift;


}


1;

__END__
