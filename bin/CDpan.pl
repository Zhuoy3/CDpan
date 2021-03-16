#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-03

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../tools/Prel/lib/";

use Getopt::Std;
use File::Spec::Functions;
use Mnet::Tee (); # This module has not modified 'say', must use 'print'
use Config::IniFiles;
use CDpan::GetPar;
use CDpan::Check;
#use CDpan::QualityControl;

my $VERSION = 'v0.0.1';
my $VERSION_TIME = 'Mar 9 2021';
my $file_par_path; #Define in advance to ensure END block can be executed

unless (@ARGV) {
    print "CDpan $VERSION ($VERSION_TIME).\n";
    print "Use option \'-h/--help\' for usage information.\n";
    exit 0;
}

# h:help; v:version; d:debug mode
our ($opt_h, $opt_v, $opt_d);
Getopt::Std::getopts('hvd') or die "Use option \'-h/--help\' for usage information.\n";

if ($opt_h) {
    HELP_MESSAGE();
    exit 0;
}

if ($opt_v) {
    print "Version $VERSION\n";
    exit 0;
}

print "Start running...\n";
print "Run in DEBUG mode.\n" if $opt_d;
print "\n====================\n\n";

# Import path of parameter file
if (@ARGV == 1) {
    ($file_par_path) = @ARGV;
}elsif (@ARGV == 0) {
    print "Please imput name of parameter file:\n";
    chomp($file_par_path = <STDIN>);
}elsif (@ARGV > 1) {
    warn "WARNING: More than two parameter files, $ARGV[1] (and any subsequent parameter files) was ignored.\n";
    ($file_par_path) = @ARGV;
}

die "ERROR: There is no such parameter file: $file_par_path.\n" unless (-e $file_par_path);
print  "\n====================\n\n";
print "Read parameters from \'$file_par_path\'.\n";

my $par = CDpan::GetPar::getpar($file_par_path);

# Check if the output folder exists
if ( -e $par->val('DATA', 'output') ) {
    my $folder_output = $par->val('DATA', 'output');
    my $folder_output_old = ( $folder_output =~ s/\/$//r ). '.old/';
    $folder_output_old =~ s/\.old\/?$/\.$$\.old\// if (-e $folder_output_old);
    rename $folder_output => $folder_output_old or die "ERROR: Cannot change the name of folder '$folder_output': $!\n";
    warn "WARNING: Folder '$folder_output' exists and has been renamed to '$old_output'";
}
mkdir $par->val('DATA', 'output') or die "ERROR: Cannot create output directory: $!\n";

# Create the folder for process file
our $folder_process = "tmp_$$/";
mkdir $folder_process or die "ERROR: Cannot create process folder '$folder_process': $!\n";

my $cwd = curdir();
chdir $folder_process or die "ERROR: Cannot chdir to '$folder_process: $!\n";

print  "\n====================\n\n";
print "Start quality control...\n";
#CDpan::QualityControl::();

print "END OF PROGRAMME.\n";
exit 0;



# HELP_MESSAGE can called by Getopt::Std when supplying '--help' option
sub HELP_MESSAGE {
    system "pod2text $FindBin::Bin/../doc/help.pod";
    exit 0;
}

END {
    my $file_log_path = $file_par_path // 'CDpan';
    $file_log_path =~ s/\.[^\.]+$//;
    Mnet::Tee::file("$file_log_path.log");
}
