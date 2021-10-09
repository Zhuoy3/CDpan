#!/usr/bin/perl

# Description:
# Author: zhuoy
# Date: 2021-10-09

package CDpanCheck;

use strict;
use warnings;

use File::Spec::Functions  qw /:ALL/;
use Config::IniFiles;
use CDpanPrint qw / :ALL /;

require Exporter;
our @ISA = qw \ Exporter \;
our @EXPORT = qw \ PreCheck \;
our %EXPORT_TAGS = ( ALL => [ @EXPORT ] );

sub PreCheck {
    (my $par) = @_;

    __CheckTools__($par);

    return 1;
}

sub __CheckTools__ {
    (my $par) = @_;

    my $exit_tools = 0;

    my @tools_needed = qw \ trim_galore
                            cutadapt
                            fastqc
                            bwa
                            gatk
                            samtools
                            masurca
                            centrifuge
                            centrifuge-kreport
                            mmseqs
                            nucmer
                            show-coords
                            RepeatMasker
                            bowtie2
                            bedtools
                            bowtie2-build \;
    print STDERR "Tools needed: @tools_needed\n";
    print STDERR "\n";

    my @tools_specify = $par->Parameters('TOOLS');
    print STDERR "Tools specified: @tools_specify\n";
    print STDERR "\n";

    print STDERR "Specified:\n";
    print STDERR "\n";
    foreach my $tools (@tools_specify) {
        unless (grep { $_ eq $tools } @tools_needed) {
            PrintWarnMessage("$tools is not needed, ignore it\n");
            $par->delval('TOOLS', $tools);
            next;
        }

        my $tools_path = $par->val('TOOLS', $tools);
        $tools_path = rel2abs($tools_path) unless file_name_is_absolute($tools_path);
        printf STDERR "%-22s", "$tools:";
        print  STDERR "$tools_path\n";

        if ( -e -x $tools_path ) {
            $par->setval('TOOLS', $tools, $tools_path);
        }
        else{
            PrintWarnMessage("$tools does not exist or lacks execute permission, search it from PATH\n");
            $par->delval('TOOLS', $tools);
        }
    }

    print STDERR "\n";
    print STDERR "Search from PATH:\n";
    print STDERR "\n";
    foreach my $tools (@tools_needed) {
        unless ( defined $par->val('TOOLS', $tools) ) {
            my $cmd_findtool = "command -v ${tools}";
            my $new_tools_path =  `$cmd_findtool`;
            chomp($new_tools_path);

            if ($new_tools_path eq '') {
                PrintErrorMessage("Couln't find tool $tools from PATH\n", 0);
                $exit_tools = 1;
            }
            else{
                printf STDERR "%-22s", "$tools:";
                print  STDERR "$new_tools_path\n";
                $par->newval('TOOLS', $tools, $new_tools_path);

            }
        }
    }

    exit(255) if $exit_tools;

    return 1;
}

1;

__END__
