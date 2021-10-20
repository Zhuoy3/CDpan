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
use CDpan::Module::Soot;
use CDpan::Module::Merge;
# use CDpan::Module::Location;
use CDpan::Module::LocationSubModule;

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

    PrintStartMessage("Start Module vot");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'vot');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Vot::Vot($par, $idv_name) or PrintErrorMessage("Module vot exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "vot" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'vot');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'vot', $output_dir);

        print STDERR "Since module vot is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'vot', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module vot");

    return 1;
};

sub Soot {
    (my $par) = @_;

    PrintStartMessage("Start Module soot");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'soot');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Soot::Soot($par, $idv_name) or PrintErrorMessage("Module soot exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "soot" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'soot');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'soot', $output_dir);

        print STDERR "Since module soot is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'soot', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module soot");

    return 1;
};

sub Merge {
    (my $par) = @_;

    PrintStartMessage("Start Module merge");

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'merge');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::Merge::Merge($par, $idv_name) or PrintErrorMessage("Module merge exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "merge" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'merge');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'merge', $output_dir);

        print STDERR "Since module merge is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'merge', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module merge");

    return 1;
};

sub Location {
    (my $par) = @_;

    PrintStartMessage("Start Module location");

    print STDERR "Processing: RepeatMasker\n";
    CDpan::Module::LocationSubModule::RepeatMasker($par)
        or PrintErrorMessage("Module location exited abnormally when do RepeatMasker");
    print STDERR "\n";

    my $work_dir = catdir($par->val('CDPAN', 'work_dir'), 'pre_location');
    mkdir $work_dir or PrintErrorMessage("Cannot create work direction $work_dir: $!");

    my $extract_dir = $par->val('RESULT', 'extract') // $par->val('LOCATION', 'extract_dir');
    $par->setval('CDPAN', 'input_dir', $extract_dir);
    my @input_idvs = sort ( File::Slurp::read_dir( $par->val('CDPAN', 'input_dir')) );

    foreach my $idv_name (@input_idvs) {
        print STDERR "Processing: $idv_name\n";
        CDpan::Module::LocationSubModule::PreLocation($par, $idv_name)
            or PrintErrorMessage("Module location exited abnormally when do RepeatMasker");

        # CDpan::Module::Location::Location($par, $idv_name) or PrintErrorMessage("Module location exited abnormally for $idv_name");
        print STDERR "\n";
    }

    if ($main::modules{ "location" }){
        my $output_dir = catdir($par->val('CDPAN', 'output_dir'), 'pre_location');
        move $work_dir, $output_dir or PrintErrorMessage("Couln't move $work_dir to $output_dir: $!");
        $par->newval('RESULT', 'location', $output_dir);

        print STDERR "Since module location is being used, program will end\n";
    }
    elsif ($main::modules{ "RUN-ALL" } or $main::modules{ "RUN-DISPLACE" }){
        $par->newval('RESULT', 'location', $work_dir);
        $par->setval('CDPAN', 'input_dir', $work_dir);

        print STDERR "Since module $main::module is being used, continue to run module align\n";
    }

    PrintEndMessage("Finish Module location");

    return 1;
};

sub RunAll;
sub RunDisplace;


#     CDpan::Integration::integration($par, $idv_name, $idv_output_folder, $folder_process)
#         or PrintErrorMessage("Operation Integration is abnormal.\n");
#     #location
#     system "python3 $FindBin::Bin/ex.py"








1;

__END__
