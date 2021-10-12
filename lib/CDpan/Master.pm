#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-10-09

package CDpan::Master;

use strict;
use warnings;

use File::Spec::Functions qw /:ALL/;
use File::Copy qw / move /;
use Config::IniFiles;
use File::Slurp;

use CDpan::Print qw / :PRINT /;

use CDpan::Module::Filter;
use CDpan::Module::Align;

# use CDpan::Extract;
# use CDpan::Assembly;
# use CDpan::Test;
# use CDpan::Judge;
# use CDpan::MMSeqs;
# use CDpan::Nucmer;
# use CDpan::DeRepeat;
# use CDpan::Recode;
# use CDpan::RepeatMasker;
# use CDpan::Align;
# use CDpan::Change;
# use CDpan::Integration;

sub Filter {
    (my $par) = @_;

    PrintStartMessage("Start Module filter");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'filter');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Filter::Filter($par, $idv_name) or PrintErrorMessage("Module filter exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "filter" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'filter');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'filter', $output_dir);

        print STDERR "Since module filter is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module filter");

    return 1;
};

sub Align {
    (my $par) = @_;

    PrintStartMessage("Start Module align");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'align');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    print STDERR "Index\n";
    CDpan::Module::Align::AlignIndex($par) or PrintErrorMessage("Module align exited abnormally when build index");
    print STDERR "\n";

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Align::Align($par, $idv_name) or PrintErrorMessage("Module align exited abnormally for $idv_name");
        print STDERR "\n";
    }

    CDpan::Module::Align::AlignRemoveIndex($par) or PrintErrorMessage("Module align exited abnormally when remove index");

    if ($main::modules{ "align" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'align');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'align', $output_dir);

        print STDERR "Since module align is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module align");

    return 1;
};

sub Extract;
sub Assembly;
sub Mope;
sub Vot;
sub Soot;
sub Merge;
sub Location;
sub RunAll;
sub RunDisplace;



# our $ref_index = 0; # mark the existence of the reference genome index used by bwa
# our $ref_dict = 0; # mark the existence of the reference genome dict used by gatkmy



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







1;

__END__
