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
use File::Copy qw / copy move /;
use File::Path qw / rmtree /;
use File::Slurp;
use File::Spec::Functions  qw /:ALL/;
use Mnet::Tee (); # This module has not modified 'say', must use 'print' for log
use Config::IniFiles;
use CDpan::Parameter;
use CDpan::QualityControl;
use CDpan::Comparison;
use CDpan::Extract;
use CDpan::Assembly;
use CDpan::Test;
use CDpan::Judge;
use CDpan::MMSeqs;
use CDpan::Nucmer;
use CDpan::DeRepeat;
use CDpan::Recode;
use CDpan::RepeatMasker;
use CDpan::Align;
use CDpan::Change;
use CDpan::Integration;

my $file_par_path = $ENV{'CDPAN_SCRIPT'} or die "Error: Cannot find parameter of script file.\n";
our $debug = $ENV{'CDPAN_DEBUG'} ? 1 : 0;
our $quiet = $ENV{'CDPAN_QUITE'} ? 1 : 0;

# couldn't work with Ment::Tee
# if ($quiet) {
#     open STDOUT, ">", "/dev/null" or die "Error: cannot redirect standard output: $!\n";
#     open STDERR, ">", "/dev/null" or die "Error: cannot redirect standard Error: $!\n";
# }

my $cwd = Cwd::getcwd;

print "CDpan v0.1.12 Jue 20 2021\n";
print "Start running...\n";
print "Debug: Run in debug mode.\n" if $debug;
print "\n================================================================================\n\n";

print "Debug: Find script file: \'$file_par_path\'.\n" if $debug;
$file_par_path = rel2abs($file_par_path);
print "Debug: Real path of script file is \'$file_par_path\'.\n" if $debug;

die "Error: There is no such parameter file: \'$file_par_path\'.\n" unless (-e $file_par_path);

print "Read parameters from \'$file_par_path\'.\n";

my $par = CDpan::Parameter::InputPar($file_par_path);

# Create the folder for process file
# our $folder_process = catdir($cwd, "tmp_$$");
# TODO 便于调试简化临时文件夹
our $folder_process = catdir($cwd, "tmp");
# if (-e $folder_process){
#     system "rm -rf $folder_process";
# }
# mkdir $folder_process or die "Error: Cannot create process folder '$folder_process': $!\n";
chdir $folder_process or die "Error: Cannot chdir to '$folder_process: $!\n";

print "\n================================================================================\n\n";
print "Start quality control...\n";

my @input_folder = sort ( File::Slurp::read_dir($par->val('DATA', 'input'), prefix => 1) );

our $ref_index = 0; # mark the existence of the reference genome index used by bwa
our $ref_dict = 0; # mark the existence of the reference genome dict used by gatkmy

my @idv_names;
foreach my $idv_folder (@input_folder) {
    my $idv_name = pop [ splitdir($idv_folder) ];
    push @idv_names, $idv_name;
    my $idv_output_folder = catdir($folder_process, $idv_name);

    #TODO　Filteration
    CDpan::QualityControl::QualityControl($par, $idv_folder, $idv_name, $idv_output_folder)
        or die "Error: Operation QualityControl is abnormal.\n";
    CDpan::Comparison::comparison($par, $idv_name, $idv_output_folder)#TODO alignment
        or die "Error: Operation Comparison is abnormal.\n";
    CDpan::Extract::extract($par, $idv_name, $idv_output_folder)
        or die "Error: Operation Extract is abnormal.\n";
    CDpan::Assembly::assembly($par, $idv_name, $idv_output_folder)
        or die "Error: Operation Assembly is abnormal.\n";
    CDpan::Test::test($par, $idv_name, $idv_output_folder)#remove contaminents
        or die "Error: Operation Test is abnormal.\n";
    CDpan::Judge::judge($par, $idv_name, $idv_output_folder)
        or die "Error: Operation Judge is abnormal.\n";
    CDpan::MMSeqs::mmseqs($par, $idv_name, $idv_output_folder)#remove redundant
        or die "Error: Operation MMSeqs is abnormal.\n";
    CDpan::Nucmer::nucmer($par, $idv_name, $idv_output_folder)#precise delete
        or die "Error: Operation Nucmer is abnormal.\n";
    CDpan::DeRepeat::de_repeat($par, $idv_name, $idv_output_folder)
        or die "Error: Operation DeRepeat is abnormal.\n";

    move "$idv_output_folder/$idv_name.filtered.mmseqs.final.fa", "$folder_process/$idv_name.fasta"
        or die "Error:Couln't move $idv_output_folder/$idv_name.filtered.mmseqs.final.fa to $folder_process/$idv_name.fasta: $!.\n";

}

CDpan::Recode::recode($folder_process, \@idv_names)
    or die "Error: Operation Recode is abnormal.\n";
CDpan::RepeatMasker::repeat_masker($par, $folder_process)#mask
    or die "Error: Operation RepeatMasker is abnormal.\n";

foreach my $idv_name (@idv_names) {
    print "\n================================================================================\n\n";
    my $idv_output_folder = catdir($folder_process, $idv_name);

    CDpan::Align::align($par, $idv_name, $idv_output_folder, $folder_process)
        or die "Error: Operation Align is abnormal.\n";
    CDpan::Change::change($par, $idv_name, $idv_output_folder)
        or die "Error: Operation Change is abnormal.\n";
    CDpan::Integration::integration($par, $idv_name, $idv_output_folder, $folder_process)
        or die "Error: Operation Integration is abnormal.\n";
    #location
    system "python3 $FindBin::Bin/ex.py"
}



# chdir $cwd;
# system "rm -rf $folder_process";

    print "\n================================================================================\n\n";
    print "END OF PROGRAMME.\n";
exit 0;

END {
    print "\n================================================================================\n\n";
    my $file_log_path = ($file_par_path // 'CDpan') =~ s/\.[^\.]+$//r;
    print "Output log in file \'$file_log_path.log\'.\n";
    Mnet::Tee::file("$file_log_path.log");
}
