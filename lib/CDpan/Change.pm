#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-07-21

package CDpan::Change;

use strict;
use warnings;
use Config::IniFiles;

sub change {
    # &extract($opt, $idv_folder_name, $output_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir) = @_;

    # Read the software path and set it to the default value
    my $samtools = $par->val('TOOLS', 'samtools');

    my $cmd_samtools1 = "$samtools view -h " .
                        "-F 256 $output_dir/$idv_folder_name.readContigAlignment.final.sam " .
                        "| $samtools sort - -n -O bam " .
                        "| $bedtools bamtobed -i stdin  " .
                        "| awk \'{OFS=\"\t\"}{print \$4,\$1,\$6,\$2,\$3}\'  " .
                        "| sort > $output_dir/$idv_folder_name.readContigAlignment.txt";
    print "Start use cmd: \'$cmd_samtools1\'.\n";
    system $cmd_samtools1
        and die "Error: Command \'$cmd_samtools1\' failed to run normally: $?\n";

    my $cmd_samtools2 = "$samtools view " .
                        "-H Sus_11_1_Links.bam " . #TODO
                        "| cat - <(awk \'FNR==NR{main[\$1]=\$0;next} \$1 in main {print main[\$1]}\'" .
                                "<($samtools view Sus_11_1_Links) " .
                                "$output_dir/$idv_folder_name.readContigAlignment.txt)" .
                        "| $samtools sort -n -O bam " .
                        "| $bedtools bamtobed -i stdin " .
                        "| awk \'{OFS=\"\t\"}{print \$4,\$1,\$6,\$2,\$3}\' " .
                        "| sed -e \'s/\/[1-2]//g\' " .
                        "| sort > $output_dir/$idv_folder_name.matchedMates.txt";
    print "Start use cmd: \'$cmd_samtools2\'.\n";
    system $cmd_samtools2
        and die "Error: Command \'$cmd_samtools2\' failed to run normally: $?\n";

    my $cmd_join = "join " .
                   "-j 1 " .
                   "$output_dir/$idv_folder_name.readContigAlignment.txt " .
                   "$output_dir/$idv_folder_name.matchedMates.txt " .
                   "> $output_dir/$idv_folder_name.mateLinks.txt";
    print "Start use cmd: \'$cmd_join\'.\n";
    system $cmd_join
        and die "Error: Command \'$cmd_join\' failed to run normally: $?\n";

}

1;

__END__
