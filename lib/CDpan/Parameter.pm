#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-07
# called by: CDpan.pl

package CDpan::Parameter;

use strict;
use warnings;
use Config::IniFiles;
use File::Spec::Functions qw /:ALL/;

sub InputPar {
    # &InputPar($opt)
    # Import parameter form file which be given by command line
    my $file_par_path = shift;

    # Since "FindBin" uses a "BEGIN" block, it'll be executed only once, and only the first caller will get it right
    # $FindBin::Bin is path to bin folder from where 'CDpan.pl' was invoked
    my $par_default = new Config::IniFiles(-file => "$FindBin::Bin/../config/par_default.ini");
    my $par = new Config::IniFiles(-file                => $file_par_path,
                                   -import              => $par_default,
                                   -allowedcommentchars => '#')
        or die "Error: Could not import parameter from \'$file_par_path\': @Config::IniFiles::errors.\n";


    #TODO 暂时调整
    #_CheckRedundant($par);
    _CheckMissing($par);
    _CheckTools($par) unless $main::debug;
    _CheckFile($par);

    #$par->WriteConfig("$file_par_path.import") if $main::debug;
    $par->WriteConfig("$file_par_path.import");

    return $par;
}

sub _GetStdPar {
    # &_GetStdPar($file)
    # Import parameter form file in the folder '/CDpan/config/'
    my $file_path = shift;
    my %output;

    open PARFILE, "<", $file_path or die "Error: Cannot open file \'$file_path\': $!";
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
                    warn "Warning: Unknown parameter found: '[$section] => $param', will be removed.\n";
                    $par->delval($section, $param);
                    $redundant = 1;
                }
            }
        }else{
            warn "Warning: Unknown parameter found: '[$section] => ...', will be removed.\n";
            $par->DeleteSection($section);
            $redundant = 1;
        }
    }

    die "Error: More unknown parameter, please check the parameter file.\n" if $redundant;
    print "Debug: Redundant check has been performed.\n" if $main::debug;

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
                print STDERR "Error: Following parameters should have values:\n" unless $findmiss;
                $findmiss = 1;
                print STDERR "Error:     $key => $value\n";
            }
        }

    }
    die "Error: Please modify the parameter file and add the expected value" if $findmiss;
    print "Debug: Missing check has been performed.\n" if $main::debug;

    return $findmiss;
}

sub _CheckTools {
    # &_CheckTools($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check whether the tools specified in par is available, and search tools in PATH if not specified
    my $par = shift;

    my @tools_needed = $par->Parameters('TOOLS');
    my @tools_missing;

    foreach my $tools (@tools_needed) {
        print "Check for $tools.\n";
        if ( defined $par->val('TOOLS', $tools) ) {
            if ( -e -x $par->val('TOOLS', $tools) ) {
                print "The $tools assign by script file seems to work fine.\n";
                _CheckVersion($tools, $par->val('TOOLS', $tools));
                next;
            }else {
                print "The $tools does not exist or lacks execute permission, delete it.\n";
                $par->delval('TOOLS', $tools);
            }
        } else {
            print "The $tools not specified.\n";
        }

        my $cmd_findtool = "command -v ${tools}";
        my $new_tools_path =  `$cmd_findtool`;

        unless ($new_tools_path =~ m/not found/) {
            chomp($new_tools_path);
            print "The $tools was found at $new_tools_path.\n";
            $par->newval('TOOLS', $tools, $new_tools_path);
            _CheckVersion($tools, $new_tools_path);
            next;
        }

        push @tools_missing, $tools;
    }

    #TODO 暂时关闭
    # die "Error: can't locate '@tools_missing', please add these to PATH or use 'TOOLS' section of parameters file.\n" .
    #     "If the above is indeed executed, please confirm whether you have execution authority.\n"
    #         if (@tools_missing);
    print "Debug: Tools check has been performed.\n" if $main::debug;
    return @tools_missing;
}

sub _CheckVersion {
    # &_CheckVersion($tools, $tool_path)
    # Output calling software information

    my $tools = shift;
    my $tool_path = shift;

    my $res;
    if ( $tools eq "trim_galore") {
        ( $res ) = `trim_galore --version` =~ m/version (\d+\.\d+.\d+)\s+/;
        print "trim_galore version: $res\n";
    }
    elsif ( $tools eq "cutadapt") {
        ( $res ) = `cutadapt --version` =~ m/(\d+\.\d+)\s+/;
        print "cutadapt version: $res\n";
    }
    elsif ( $tools eq "fastqc") {
        ( $res ) = `fastqc --version` =~ m/FastQC v(\d+\.\d+.\d+)\s+/;
        print "fastqc version: $res\n";
    }
    elsif ( $tools eq "bwa") {
    #     ( $res ) = `bwa` =~ m/Version: (\d+\.\d+.\d+)-/;
    #    print "bwa version: $res\n";
    # this check for bwa seem couldn't work normal
    }
    elsif ( $tools eq "gatk") {
        # open STDERR, ">", "/dev/null" or die "Error: cannot redirect standard Error: $!";
        # ( $res ) = `gatk --version` =~ m/\(GATK\) v(\d+\.\d+.\d+.\d+)\s+/;
        # print "gatk version: $res\n";
        # open STDERR, ">&$main::cp_stderr" or die "Error: Can't restore stderr: $!";
    }
    elsif ( $tools eq "samtools") {
        ( $res ) = `samtools --version` =~ m/samtools (\d+\.\d+)\s+/;
        print "samtools version: $res\n";
    }
    elsif ( $tools eq "masurca") {
        # ( $res ) = `masurca --version` =~ m/version (\d+\.\d+.\d+)\s+/;
        # print "masurca version: $res\n";
        # this check for masurca seem couldn't work normal
    }
    elsif ( $tools eq "mmseqs") {
        ;
        # this check for masurca seem couldn't work normal
    }
    elsif ( $tools eq "RepeatMasker") {
        ;
        # this check for masurca seem couldn't work normal
    }
    else {
        warn "Warning: Unknown software: $tools\n";
    }

    print "Debug: $tools is $tool_path\n" if $main::debug;

    return 1;
}

sub _CheckFile {
    # &checkpar_tool($opt)
    # $opt is a quotation in 'Config::IniFiles' format
    # Check if the required file exists
    my $par = shift;

    unless ( -e -r -d $par->val('DATA', 'input') ){
        die "Error: Could not open input folder: " . $par->val('DATA', 'input') . "\n";
    }

    unless ( -e -r $par->val('DATA', 'ref') ){
        die "Error: Could not open ref file: " . $par->val('DATA', 'ref') . "\n";
    }
    #TODO 如果太多可以写成循环
    unless ( -e -r $par->val('DATA', 'qry') ){
        die "Error: Could not open ref file: " . $par->val('DATA', 'qry') . "\n";
    }

    # Check if the output folder exists
    if ( -e $par->val('DATA', 'output') ) {
        my @dir_files = File::Slurp::read_dir($par->val('DATA', 'output'), prefix => 1);
        if ( @dir_files ) {
            print "Output folder " . $par->val('DATA', 'output') . " exists and has file.\n";
            my $folder_output = $par->val('DATA', 'output');
            my $folder_output_old = ( $folder_output =~ s/\/$//r ). '.old/';
            $folder_output_old =~ s/\.old\/?$/\.$$\.old\// if (-e $folder_output_old);
            rename $folder_output => $folder_output_old or die "Error: Cannot change the name of folder '$folder_output': $!\n";
            warn "Warning: Folder '$folder_output' exists and has been renamed to '$folder_output_old'";
            mkdir $par->val('DATA', 'output') or die "Error: Cannot create output folder: $!\n";
            print "Output folder " . $par->val('DATA', 'output') . " exists was created.\n";
        } else {
            print "Output folder " . $par->val('DATA', 'output') . " exists but is empty.\n";
        }
    } else {
        mkdir $par->val('DATA', 'output') or die "Error: Cannot create output folder: $!\n";
        print "Output folder " . $par->val('DATA', 'output') . " was created.\n";
    }


    print "Debug: Files check has been performed.\n" if $main::debug;

    return 1;
}


1;

__END__
