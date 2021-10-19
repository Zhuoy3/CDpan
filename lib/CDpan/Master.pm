#!/usr/bin/env perl

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
use CDpan::Module::Extract;
use CDpan::Module::Assembly;
use CDpan::Module::Mope;
use CDpan::Module::Vot;

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
        $par->newval('RESULT', 'filter', $work_dir);
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
        $par->newval('RESULT', 'align', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module align");

    return 1;
};

sub Extract {
    (my $par) = @_;

    PrintStartMessage("Start Module extract");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'extract');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Extract::Extract($par, $idv_name) or PrintErrorMessage("Module extract exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "extract" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'extract');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'extract', $output_dir);

        print STDERR "Since module extract is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'extract', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module extract");

    return 1;
};

sub Assembly {
    (my $par) = @_;

    PrintStartMessage("Start Module assembly");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'assembly');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Assembly::Assembly($par, $idv_name) or PrintErrorMessage("Module assembly exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "assembly" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'assembly');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'assembly', $output_dir);

        print STDERR "Since module assembly is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'assembly', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module assembly");

    return 1;
};

sub Mope {
    (my $par) = @_;

    PrintStartMessage("Start Module mope");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'mope');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Mope::Mope($par, $idv_name) or PrintErrorMessage("Module mope exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "mope" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'mope');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'mope', $output_dir);

        print STDERR "Since module mope is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'mope', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module mope");

    return 1;
};


sub Vot {
    (my $par) = @_;

    PrintStartMessage("Start Module mope");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'mope');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Mope::Mope($par, $idv_name) or PrintErrorMessage("Module mope exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "mope" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'mope');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'mope', $output_dir);

        print STDERR "Since module mope is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'mope', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module mope");

    return 1;
};

sub Soot;
sub Merge;
sub Location;
sub RunAll;
sub RunDisplace;


#     CDpan::MMSeqs::mmseqs($par, $idv_name, $idv_output_folder)#remove redundant
#         or PrintErrorMessage("Operation MMSeqs is abnormal.\n");
#     CDpan::Nucmer::nucmer($par, $idv_name, $idv_output_folder)#precise delete
#         or PrintErrorMessage("Operation Nucmer is abnormal.\n");
#     CDpan::DeRepeat::de_repeat($par, $idv_name, $idv_output_folder)
#         or PrintErrorMessage("Operation DeRepeat is abnormal.\n");

#     move "$idv_output_folder/$idv_name.filtered.mmseqs.final.fa", "$folder_process/$idv_name.fasta"
#         or PrintErrorMessage("Couln't move $idv_output_folder/$idv_name.filtered.mmseqs.final.fa to $folder_process/$idv_name.fasta: $!.\n");


# CDpan::Recode::recode($folder_process, \@idv_names)
#     or PrintErrorMessage("Operation Recode is abnormal.\n");
# CDpan::RepeatMasker::repeat_masker($par, $folder_process)#mask
#     or PrintErrorMessage("Operation RepeatMasker is abnormal.\n");

# foreach my $idv_name (@idv_names) {
#     print "\n================================================================================\n\n";
#     my $idv_output_folder = catdir($folder_process, $idv_name);

#     CDpan::Align::align($par, $idv_name, $idv_output_folder, $folder_process)
#         or PrintErrorMessage("Operation Align is abnormal.\n");
#     CDpan::Change::change($par, $idv_name, $idv_output_folder)
#         or PrintErrorMessage("Operation Change is abnormal.\n");
#     CDpan::Integration::integration($par, $idv_name, $idv_output_folder, $folder_process)
#         or PrintErrorMessage("Operation Integration is abnormal.\n");
#     #location
#     system "python3 $FindBin::Bin/ex.py"
# }







1;

__END__
