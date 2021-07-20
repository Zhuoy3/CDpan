#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-03-03

use strict;
use warnings;
use feature qw /say/;

use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../tools/Perl/lib/";

use Cwd;
use File::Slurp;
use File::Spec::Functions  qw /:ALL/;
use Mnet::Tee (); # This module has not modified 'say', must use 'print' for log
use Config::IniFiles;
use CDpan::Parameter;
use CDpan::QualityControl;
use CDpan::Comparison;
use CDpan::Extract;
use CDpan::Assembly;

my $file_par_path = $ENV{'CDPAN_SCRIPT'} or die "ERROR: Cannot find parameter of script file.\n";
our $debug = $ENV{'CDPAN_DEBUG'} ? 1 : 0;
our $quiet = $ENV{'CDPAN_QUITE'} ? 1 : 0;

if ($quiet) {
    open STDOUT, ">", "/dev/null" or die "ERROR: cannot redirect standard output";
    open STDERR, ">", "/dev/null" or die "ERROR: cannot redirect standard error";
}

my $cwd = Cwd::getcwd;

print "CDpan v0.1.12 Jue 20 2021\n";
print "Start running...\n";
print "DEBUG: Run in DEBUG mode.\n" if $debug;
print "\n================================================================================\n\n";

print "DEBUG: Find script file: \'$file_par_path\'.\n" if $debug;
$file_par_path = rel2abs($file_par_path);
print "DEBUG: Real path of script file is \'$file_par_path\'.\n" if $debug;

die "ERROR: There is no such parameter file: \'$file_par_path\'.\n" unless (-e $file_par_path);

print "Read parameters from \'$file_par_path\'.\n";

my $par = CDpan::Parameter::InputPar($file_par_path);

# # Check if the output folder exists
# if ( -e $par->val('DATA', 'output') ) {
#     my $folder_output = $par->val('DATA', 'output');
#     my $folder_output_old = ( $folder_output =~ s/\/$//r ). '.old/';
#     $folder_output_old =~ s/\.old\/?$/\.$$\.old\// if (-e $folder_output_old);
#     rename $folder_output => $folder_output_old or die "ERROR: Cannot change the name of folder '$folder_output': $!\n";
#     warn "WARNING: Folder '$folder_output' exists and has been renamed to '$folder_output_old'";
# }
# mkdir $par->val('DATA', 'output') or die "ERROR: Cannot create output directory: $!\n";

# # Create the folder for process file
# our $folder_process = catdir($cwd, "tmp_$$");
# mkdir $folder_process or die "ERROR: Cannot create process folder '$folder_process': $!\n";
# chdir $folder_process or die "ERROR: Cannot chdir to '$folder_process: $!\n";

# chdir $cwd;
# #system "rm -rf $folder_process";

# print  "\n====================\n\n";
# print "Start quality control...\n";

# #TODO 这一块可以使用子进程实现多线程？ pro.pl

# my $input_directory = $par->val('DATA', 'input');
# die "ERROR: Data folder \'$input_directory\' does not exist." unless ( -e $input_directory );
# my @input_folder = sort ( File::Slurp::read_dir($input_directory, prefix => 1) );

# our $ref_index = 0; # mark the existence of the reference genome index used by bwa
# our $ref_dict = 0; # mark the existence of the reference genome dict used by gatkmy

# foreach my $idv_folder (@input_folder) {
#     my $idv_name = pop [ splitdir($idv_folder) ];
#     my $idv_output_folder = catdir($folder_process, $idv_name);

#     CDpan::QualityControl::QualityControl($par, $idv_folder, $idv_name, $idv_output_folder)
#         or die "ERROR: Operation QualityControl is abnormal.\n";
#     CDpan::Comparison::comparison($par, $idv_name, $idv_output_folder)
#         or die "ERROR: Operation Comparison is abnormal.\n";
#     CDpan::Extract::extract($par, $idv_name, $idv_output_folder)
#         or die "ERROR: Operation Extract is abnormal.\n";
#     CDpan::Assembly::assembly($par, $idv_name, $idv_output_folder)
#         or die "ERROR: Operation Assembly is abnormal.\n";
# }

print "\n================================================================================\n\n";
print "END OF PROGRAMME.\n";
exit 0;

END {
    my $file_log_path = ($file_par_path // 'CDpan') =~ s/\.[^\.]+$//r;
    print "Output log in file \'$file_log_path.log\'.\n";
    Mnet::Tee::file("$file_log_path.log");
}
