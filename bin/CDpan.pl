#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-03

use strict;
use warnings;

use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../tools/Perl/lib/";

use Cwd;
use Getopt::Std;
use File::Slurp;
use File::Spec::Functions  qw /:ALL/;
use Mnet::Tee (); # This module has not modified 'say', must use 'print'
use Config::IniFiles;
use CDpan::GetPar;
use CDpan::Check;
#use CDpan::GetSample;
use CDpan::QualityControl;
use CDpan::Comparison;
use CDpan::Extract;
use CDpan::Assembly;

my $VERSION      = 'v0.1.8';
my $VERSION_TIME = 'May 27 2021';
my $file_par_path; #Define in advance to ensure END block can be executed

unless (@ARGV) {
    print "CDpan $VERSION ($VERSION_TIME).\n";
    print "Use option \'-h/--help\' for usage information.\n";
    exit 0;
}
#TODO 参数解析不行换用python或是C++
# h:help; v:version; d:debug mode
our ($opt_h, $opt_v, $opt_d);
Getopt::Std::getopts('hvd') or die "Error parameters \'@ARGV\', please use option \'-h/--help\' for usage information.\n";

if ($opt_h) {
    HELP_MESSAGE();
    exit 0;
}

if ($opt_v) {
    print "Version $VERSION\n";
    exit 0;
}

my $cwd = Cwd::getcwd;
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

$file_par_path = rel2abs($file_par_path);
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
    warn "WARNING: Folder '$folder_output' exists and has been renamed to '$folder_output_old'";
}
mkdir $par->val('DATA', 'output') or die "ERROR: Cannot create output directory: $!\n";

# Create the folder for process file
our $folder_process = catdir($cwd, "tmp_$$");
mkdir $folder_process or die "ERROR: Cannot create process folder '$folder_process': $!\n";
chdir $folder_process or die "ERROR: Cannot chdir to '$folder_process: $!\n";

chdir $cwd;
#system "rm -rf $folder_process";

print  "\n====================\n\n";
print "Start quality control...\n";

#TODO 这一块可以使用子进程实现多线程？ pro.pl

my $input_directory = $par->val('DATA', 'input');
die "ERROR: Data folder \'$input_directory\' does not exist." unless ( -e $input_directory );
my @input_folder = sort ( File::Slurp::read_dir($input_directory, prefix => 1) );

our $ref_index = 0; # mark the existence of the reference genome index used by bwa
our $ref_dict = 0; # mark the existence of the reference genome dict used by gatkmy

foreach my $idv_folder (@input_folder) {
    my $idv_name = pop [ splitdir($idv_folder) ];
    my $idv_output_folder = catdir($folder_process, $idv_name);

    CDpan::QualityControl::QualityControl($par, $idv_folder, $idv_name, $idv_output_folder)
        or die "ERROR: Operation QualityControl is abnormal.\n";
    CDpan::Comparison::comparison($par, $idv_name, $idv_output_folder)
        or die "ERROR: Operation Comparison is abnormal.\n";
    CDpan::Extract::extract($par, $idv_name, $idv_output_folder)
        or die "ERROR: Operation Extract is abnormal.\n";
    CDpan::Assembly::assembly($par, $idv_name, $idv_output_folder)
        or die "ERROR: Operation Assembly is abnormal.\n";
}

#TODO
print "END OF PROGRAMME.\n";
exit 0;

# HELP_MESSAGE can called by Getopt::Std when supplying '--help' option
sub HELP_MESSAGE {
    system "pod2text $FindBin::Bin/../doc/help.pod";
    exit 0;
}

END {
    unless ($opt_h or $opt_v) {
        my $file_log_path = $file_par_path // 'CDpan';
        $file_log_path =~ s/\.[^\.]+$//;
        print "Output log in file \'$file_log_path.log\'.\n";
        Mnet::Tee::file("$file_log_path.log");
    }
}
