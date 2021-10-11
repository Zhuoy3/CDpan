#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-03-03

use strict;
use warnings;
use feature qw\ say \;

use FindBin;
use lib "$FindBin::Bin/../lib";
use lib "$FindBin::Bin/../tools/Perl/lib/";

use Cwd;
# use File::Copy qw / copy move /;
# use File::Path qw / rmtree /;
use File::Spec::Functions  qw /:ALL/;
use Config::IniFiles;

use CDpanPrint qw / :ALL /;
use CDpanCheck qw / :ALL /;

#-------------------------------------------------------------------------------
#------------------------------------ START ------------------------------------
#-------------------------------------------------------------------------------
our $version = "0.2.2";
our $date = "Oct 8 2021";

Usage() if @ARGV == 0;

#-------------------------------------------------------------------------------
#---------------------------------- CMD LINE -----------------------------------
#-------------------------------------------------------------------------------
our $cwd = Cwd::getcwd;
our $datestring = localtime();
our $module = shift @ARGV;

if ( $module =~ m/^-/ ) {
    if ( $module eq '-v' or $module eq '--version') {
        PrintExitMessage("CDpan $version $date\n");
    }elsif ( $module eq '-h' or $module eq '--help') {
        Usage();
    }
}

our %modules = (
    "filter"            => 0,
    "align"             => 0,
    "extract"           => 0,
    "assembly"          => 0,
    "mope"              => 0,
    "judge"             => 0,
    "vot"               => 0,
    "soot"              => 0,
    "merge"             => 0,
    "RUN-ALL"           => 0,
    "RUN-DISPLACE"      => 0,
    "MAKE-CONFIG"       => 0,
);

if ( defined $modules{$module} ) {
    $modules{$module} = 1;
} else {
    PrintExitMessage("Invalid module: $module\n");
}

our $input_dir;
our $config_file;
our $output_prefix;
our $output_dir = $cwd;
our $work_dir = catdir($output_dir,'./cdpan_tmp');
our $save_process = 0;
our $no_quality_control = 0;

while (@ARGV) {
    my $option = shift;

    if    ( $option eq '-i' or $option eq '--input') {
        $input_dir = shift;
        PrintExitMessage("Command \'$option\' is missing parameters\n") if ( ($input_dir // '-' )=~ /^-/ );
        $input_dir = rel2abs($input_dir) unless file_name_is_absolute($input_dir);
    }
    elsif ( $option eq '-c' or $option eq '--config') {
        $config_file = shift;
        PrintExitMessage("Command \'$option\' is missing parameters\n") if ( ($config_file // '-' )=~ /^-/ );
        $config_file = rel2abs($config_file) unless file_name_is_absolute($config_file);
    }
    elsif ( $option eq '-w' or $option eq '--work_dir') {
        $work_dir = shift;
        PrintExitMessage("Command \'$option\' is missing parameters\n") if ( ($work_dir // '-' )=~ /^-/ );
        $work_dir = rel2abs($work_dir) unless file_name_is_absolute($work_dir);
    }
    elsif ( $option eq '-o' or $option eq '--output') {
        $output_prefix = shift;
        PrintExitMessage("Command \'$option\' is missing parameters\n") if ( ($output_prefix // '-' )=~ /^-/ );
    }
    elsif ( $option eq '-O' or $option eq '--output_dir') {
        $output_dir  = shift;
        PrintExitMessage("Command \'$option\' is missing parameters\n") if ( ($output_dir // '-' )=~ /^-/ );
        $output_dir = rel2abs($output_dir) unless file_name_is_absolute($output_dir);
    }
    elsif ( $option eq '-p' or $option eq '--process') {
        $save_process = 1;
    }
    elsif (                    $option eq '--no-qc') {
        $no_quality_control = 1;
    }
    elsif ( $option eq '-v' or $option eq '--version') {
        print STDERR "CDpan $version $date\n";
        PrintExitMessage();
    }
    elsif ( $option eq '-h' or $option eq '--help') {
        Usage();
    }
    else {
        PrintExitMessage("Invalid command: $option\n");
    }
}

unless ( defined $input_dir ) {
    PrintExitMessage("Parameter \'input_dir\' is required\n");
}
# unless ( not defined $config_file ) { #TODO and ( $modules{$module} or $modules{$module}) ) {
#     PrintExitMessage("Parameter \'config_file\' is required\n");
# }
unless ( defined $output_prefix ) {
    ( undef, undef, $output_prefix ) = splitpath($input_dir);
}

#-------------------------------------------------------------------------------
#------------------------------------ MAIN -------------------------------------
#-------------------------------------------------------------------------------

my $main_save_process       = $save_process?"True":"False";
my $main_no_quality_control = $no_quality_control?"True":"False";

print STDERR "
                     _____  ______
                    /  __ \\ |  _  \\
                    | /  \\/ | | | |  _ __     __ _   _ __
                    | |     | | | | | \'_ \\   / _` | | \'_ \\
                    | \\__/\\ | |/ /  | |_) | | (_| | | | | |
                    \\____/  |___/   | .__/   \\__,_| |_| |_|
                                    | |
                                    |_|

--------------------------------------------------------------------------------

Version:                 $version ( $date )
Runing:                  $datestring

Module :                 $module

Input directory:         $input_dir
Config file:             $config_file
Work directory:          $work_dir
Prefix of output file:   $output_prefix
Output directory:        $output_dir

Save process file:       $main_save_process
No quality control:      $main_no_quality_control

";

#-------------------------------------------------------------------------------
#----------------------------------- CHECK -------------------------------------
#-------------------------------------------------------------------------------

#TODO read input file

our $par = new Config::IniFiles(-file => $config_file, -allowedcommentchars => '#')
    or PrintErrorMessage("Could not read config file from \'$config_file\': @Config::IniFiles::errors\n");
PreCheck($par);

#TODO 暂时调整
# #_CheckRedundant($par);
# _CheckMissing($par);
# _CheckTools($par) unless $main::debug;
# _CheckFile($par);

# #$par->WriteConfig("$file_par_path.import") if $main::debug;
# $par->WriteConfig("$file_par_path.import");




# $file_par_path = rel2abs($file_par_path);



# my $par = CDpan::Parameter::InputPar($file_par_path);

# # Create the folder for process file
# # our $folder_process = catdir($cwd, "tmp_$$");
# # TODO 便于调试简化临时文件夹
# our $folder_process = catdir($cwd, "tmp");
# # if (-e $folder_process){
# #     system "rm -rf $folder_process";
# # }
# # mkdir $folder_process or die "Error: Cannot create process folder '$folder_process': $!\n";
# chdir $folder_process or die "Error: Cannot chdir to '$folder_process: $!\n";

# print "\n================================================================================\n\n";
# print "Start quality control...\n";

# my @input_folder = sort ( File::Slurp::read_dir($par->val('DATA', 'input'), prefix => 1) );

# our $ref_index = 0; # mark the existence of the reference genome index used by bwa
# our $ref_dict = 0; # mark the existence of the reference genome dict used by gatkmy

# my @idv_names;
# foreach my $idv_folder (@input_folder) {
#     my $idv_name = pop [ splitdir($idv_folder) ];
#     push @idv_names, $idv_name;
#     my $idv_output_folder = catdir($folder_process, $idv_name);

#     #TODO　Filteration
#     CDpan::QualityControl::QualityControl($par, $idv_folder, $idv_name, $idv_output_folder)
#         or die "Error: Operation QualityControl is abnormal.\n";
#     CDpan::Comparison::comparison($par, $idv_name, $idv_output_folder)#TODO alignment
#         or die "Error: Operation Comparison is abnormal.\n";
#     CDpan::Extract::extract($par, $idv_name, $idv_output_folder)
#         or die "Error: Operation Extract is abnormal.\n";
#     CDpan::Assembly::assembly($par, $idv_name, $idv_output_folder)
#         or die "Error: Operation Assembly is abnormal.\n";
#     CDpan::Test::test($par, $idv_name, $idv_output_folder)#remove contaminents
#         or die "Error: Operation Test is abnormal.\n";
#     CDpan::Judge::judge($par, $idv_name, $idv_output_folder)
#         or die "Error: Operation Judge is abnormal.\n";
#     CDpan::MMSeqs::mmseqs($par, $idv_name, $idv_output_folder)#remove redundant
#         or die "Error: Operation MMSeqs is abnormal.\n";
#     CDpan::Nucmer::nucmer($par, $idv_name, $idv_output_folder)#precise delete
#         or die "Error: Operation Nucmer is abnormal.\n";
#     CDpan::DeRepeat::de_repeat($par, $idv_name, $idv_output_folder)
#         or die "Error: Operation DeRepeat is abnormal.\n";

#     move "$idv_output_folder/$idv_name.filtered.mmseqs.final.fa", "$folder_process/$idv_name.fasta"
#         or die "Error:Couln't move $idv_output_folder/$idv_name.filtered.mmseqs.final.fa to $folder_process/$idv_name.fasta: $!.\n";

# }

# CDpan::Recode::recode($folder_process, \@idv_names)
#     or die "Error: Operation Recode is abnormal.\n";
# CDpan::RepeatMasker::repeat_masker($par, $folder_process)#mask
#     or die "Error: Operation RepeatMasker is abnormal.\n";

# foreach my $idv_name (@idv_names) {
#     print "\n================================================================================\n\n";
#     my $idv_output_folder = catdir($folder_process, $idv_name);

#     CDpan::Align::align($par, $idv_name, $idv_output_folder, $folder_process)
#         or die "Error: Operation Align is abnormal.\n";
#     CDpan::Change::change($par, $idv_name, $idv_output_folder)
#         or die "Error: Operation Change is abnormal.\n";
#     CDpan::Integration::integration($par, $idv_name, $idv_output_folder, $folder_process)
#         or die "Error: Operation Integration is abnormal.\n";
#     #location
#     system "python3 $FindBin::Bin/ex.py"
# }



# chdir $cwd;
# system "rm -rf $folder_process";

#     print "\n================================================================================\n\n";
#     print "END OF PROGRAMME.\n";
# exit 0;

sub Usage {

print STDERR "
Program: CDpan
Version: $version

Usage: CDpan <module> -i input_dir [options]

Note: To use CDpan, you need to specify a module and a directory of input file. Each sub-directory
      of the input directory should be a separate individual and contain file which will be used
      by the module. CDpan will interpret directory name as individual names. For more information,
      please refer to the manual from https://github.com/kimi-du-bio/CDpan.

Module: filter
        align
        extract
        assembly
        mpope
        judge
        vot
        soot
        merge
        location

        RUN-ALL         Run all modules of CDpan
        RUN-DISPLACE    Run all modules of CDpan except location

Options: -i, --input         path of input directory          ( Mandatory )
         -c, --config        path of config file              ( Required by )
         -w, --work_dir      path of work directory           ( Default is \'\${output_dir}/cdpan_tmp\' )
         -p, --process       whether to keep process files    ( Default is False )
         -o, --output        prefix of the output file        ( Default is name of input directory )
         -O, --output_dir    output directory                 ( Default is the current directory )
             --no-qc         no quality control               ( Default is False )
         -v, --version       print version message
         -h, --help          print help message

";
#TODO
exit(-1);
}
