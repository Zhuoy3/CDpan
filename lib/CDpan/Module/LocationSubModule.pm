#!/usr/bin/env perl

# Description:
# Author: zhuoy
# Date: 2021-08-11

package CDpan::Module::LocationSubModule;

use strict;
use warnings;

use Config::IniFiles;
use File::Spec::Functions qw /:ALL/;
use File::Copy qw / copy /;

use CDpan::Print qw / :PRINT /;

sub RepeatMasker {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    (my $par) = @_;

    my $output_dir = catdir($par->val('CDPAN', 'work_dir'), 'repeat_masker');
    mkdir $output_dir or PrintErrorMessage("Cannot create work direction $output_dir: $!");

    my $input_file = catfile($par->val('CDPAN', 'input_dir'), 'merge.fasta');
    if ($main::modules{ "location" }){
        unless ( -e $input_file){
            PrintErrorMessage("The input file $input_file does not exist, whether the input direction is the output direction of Module merge");
        }
    }

    copy $input_file, "$output_dir/merge.fasta"
        or PrintErrorMessage("Cannot copy file $input_file to $output_dir/merge.fasta: $!");

    # Read the software path
    my $repeat_masker = $par->val('TOOLS', 'RepeatMasker');

    my $thread  = $par->val('CDPAN',    'thread');
    my $species = $par->val('LOCATION', 'species');

    my $cmd_repeat_masker = "$repeat_masker " .
                            "-parallel $thread " .
                            "-nolow " .
                            "-species $species " .
                            "$output_dir/merge.fasta " .
                            "> $output_dir/repeat_masker.log";
    # print "Start use cmd: \'$cmd_repeat_masker\'.\n";
    PrintProcessMessage('screen sequence to %%*', "$output_dir/merge.fasta");
    system $cmd_repeat_masker
        and PrintErrorMessage("Command \'$cmd_repeat_masker\' failed to run normally: $?\n");

    # Read the software path
    my $samtools = $par->val('TOOLS', 'samtools');

    my $cmd_samtools = "$samtools " .
                       "faidx " .
                       "$output_dir/merge.fasta";
    # print "Start use cmd: \'$cmd_samtools\'.\n";
    PrintProcessMessage('queries regions to %%', "$output_dir/merge.fasta.fai");
    system $cmd_samtools
        and PrintErrorMessage("Command \'$cmd_samtools\' failed to run normally: $?\n");

    # Read the software path
    my $bowtie2_build = $par->val('TOOLS', 'bowtie2-build');
    my $cmd_bowtie2_build = "$bowtie2_build " .
                            "$output_dir/merge.fasta.masked " .
                            "$output_dir/index " .
                            "> $output_dir/bowtie2_build1.log " .
                            "2> $output_dir/bowtie2_build2.log";
    # print "Start use cmd: \'$cmd_bowtie2_build\'.\n";
    PrintProcessMessage('build index to %%*', "$output_dir/index");
    system $cmd_bowtie2_build
        and PrintErrorMessage("Command \'$cmd_bowtie2_build\' failed to run normally: $?\n");

    $par->newval('RESULT', 'repeat_masker', $output_dir);

    return 1;
}

sub PreLocation {
    (my $par, my $idv_name) = @_;

    my $output_dir = catdir($par->val('CDPAN', 'work_dir'),'pre_location', $idv_name);
    mkdir $output_dir or PrintErrorMessage("Cannot create direction $output_dir: $!");

    my $output_file_prefix = catfile($output_dir, $idv_name);

    my $repeat_masker_dir = $par->val('RESULT', 'repeat_masker');

    my $input_file_prefix = catfile($par->val('CDPAN', 'input_dir'), $idv_name, $idv_name);
    if ($main::modules{ "location" }){
        unless ( -e "${input_file_prefix}.singleUnmapped_R2.fq"){
            PrintErrorMessage("The input file ${input_file_prefix}.singleUnmapped_R2.fq does not exist, whether the input direction is the output direction of Module vot");
        }
    }

    # Read the software path
    my $bowtie2 = $par->val('TOOLS',  'bowtie2');
    my $thread  = $par->val('CDPAN',  'thread');

    my $cmd_align = "$bowtie2 " .
                    "-x $repeat_masker_dir/index " .
                    "-U ${input_file_prefix}.singleUnmapped_R2.fq,${input_file_prefix}.singleUnmapped_R2.fq " .
                    "-S ${output_file_prefix}.readContigAlignment.final.sam " .
                    "-p $thread " .
                    "2> ${output_file_prefix}.bowtie2.log";
    # print "Start use cmd: \'$cmd_align\'.\n";
    PrintProcessMessage('aligning sequencing to %%*', "${output_file_prefix}.readContigAlignment.final.sam");
    system $cmd_align
        and PrintErrorMessage("Command \'$cmd_align\' failed to run normally: $?\n");

    (my $par, my $idv_folder_name, my $output_dir) = @_;

    # Read the software path
    my $samtools = $par->val('TOOLS', 'samtools');
    my $bedtools = $par->val('TOOLS', 'bedtools');

    my $cmd_samtools1 = "$samtools view -@ $thread -h " .
                        "-F 256 ${output_file_prefix}.readContigAlignment.final.sam " .
                        "| $samtools sort - -n -O bam -@ $thread " .
                        "2> /dev/null " .
                        "| $bedtools bamtobed -i stdin " .
                        "| awk \'{OFS=\"\\t\"}{print \$4,\$1,\$6,\$2,\$3}\' " .
                        "| sort > ${output_file_prefix}.readContigAlignment.txt";
    # print "Start use cmd: \'$cmd_samtools1\'.\n";
    PrintProcessMessage('change to %%', "${output_file_prefix}.readContigAlignment.txt");
    system $cmd_samtools1
        and PrintErrorMessage("Command \'$cmd_samtools1\' failed to run normally: $?\n");

    my $cmd_samtools2 = "$samtools view -@ $thread " .
                        "-H ${input_file_prefix}.Sus_11_1_Links.bam " .
                        "| cat - <(awk \'FNR==NR{main[\$1]=\$0;next} \$1 in main {print main[\$1]}\' " .
                                "<($samtools view -@ $thread ${input_file_prefix}.Sus_11_1_Links.bam) " .
                                "${output_file_prefix}.readContigAlignment.txt) " .
                        "| $samtools sort -n -O bam -@ $thread " .
                        "| $bedtools bamtobed -i stdin " .
                        "| awk \'{OFS=\"\\t\"}{print \$4,\$1,\$6,\$2,\$3}\' " .
                        "| sed -e \'s/\\/[1-2]//g\' " .
                        "| sort > ${output_file_prefix}.matchedMates.txt";
    # print "Start use cmd: \'$cmd_samtools2\'.\n";
    PrintProcessMessage('change to %%', "${output_file_prefix}.matchedMates.txt");

    open my $TMP, '>', "${output_file_prefix}.tmp.sh";
    print $TMP $cmd_samtools2;
    close $TMP;

    system "bash ${output_file_prefix}.tmp.sh"
        and PrintErrorMessage("Command \'$cmd_samtools2\' failed to run normally: $?\n");

    unlink "${output_file_prefix}.tmp.sh";

    my $cmd_join = "join " .
                   "-j 1 " .
                   "${output_file_prefix}.readContigAlignment.txt " .
                   "${output_file_prefix}.matchedMates.txt " .
                   "> ${output_file_prefix}.mateLinks.txt";
    # print "Start use cmd: \'$cmd_join\'.\n";
    PrintProcessMessage('join to %%', "${output_file_prefix}.mateLinks.txt");
    system $cmd_join
        and PrintErrorMessage("Command \'$cmd_join\' failed to run normally: $?\n");

    return 1;
}

1;

__END__
