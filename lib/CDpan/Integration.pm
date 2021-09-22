#!/usr/bin/perl

# Description:
# Author: Zhuo Yue
# Date: 2021-08-21

package CDpan::Integration;

use strict;
use warnings;
use Config::IniFiles;
use File::Copy qw / copy move /;
use List::Util qw / max min /;
use File::Path qw / mkpath /;

sub integration {
    # &extract($opt, $idv_folder_name, $output_dir, $work_dir)
    # $opt is a quotation in 'Config::IniFiles' format
    # $idv_folder_name is the name of the individual
    # $output_dir is a folder path which is used to output
    (my $par, my $idv_folder_name, my $output_dir, my $work_dir) = @_;

    mkdir "$output_dir/link_new"
        or die "Error: Cannot create process folder '$output_dir/link_new': $!\n";

    # switch to destination file
    copy "$output_dir/$idv_folder_name.mateLinks.txt", "$output_dir/link_new/mateLinks.txt"
        or die "Error: Cannot copy file '$output_dir/$idv_folder_name.mateLinks.txt': $!\n";
    copy "$work_dir/all.fasta.fai", "$output_dir/link_new/all.fasta.fai"
        or die "Error: Cannot copy file '$work_dir/all.fasta.fai': $!\n";

    chdir "$output_dir/link_new"
        or die "Error: Cannot chdir to '$output_dir/link_new': $!\n";

	mkdir "./1";
	mkdir "./2";
	mkdir "./3";
	mkdir "./4";

    # archive the contig name and length
    system "awk '{print \$1,\$2}' ./all.fasta.fai > ./contig.name"
        and die "Error: Command failed to run normally: $?\n";

    # 1do
    open my $CONTIG, '<', "./mateLinks.txt"
        or die "Error: Cannot open file '$output_dir/link_new/mateLinks.txt': $!\n";

    my %contig;
    my $line = 0;
    while(<$CONTIG>) {
        chomp;
        my @contig_lines = split /\s+/, $_;
        $contig{$contig_lines[1]}{$line} = $_;
        $line++;
    }
    close $CONTIG;


    open my $INPUT, "<", "./contig.name"
        or die "Error: Cannot open file '$output_dir/link_new/contig.name': $!\n";
    open my $OUTPUT1, ">", "./1/1.name"
        or die "Error: Cannot create file '$output_dir/link_new/1/1.name': $!\n";
    open my $OUTPUT2, ">", "./1/2to4.name"
        or die "Error: Cannot create file '$output_dir/link_new/1/2to4.name': $!\n";

    my @contig_2to4;
    while(<$INPUT>) {
        chomp;
        my @lines_links = split /\s+/, $_;
        my $have_contigs = 0;

        open my $LINK, ">", "./4/$lines_links[0].link"
            or die "Error: Cannot create file '$output_dir/link_new/4/$lines_links[0].link': $!\n";

        if( exists $contig{$lines_links[0]} ) {
            foreach my $key (keys %{$contig{$lines_links[0]}}) {
                print $LINK "$contig{$lines_links[0]}{$key}\n";
            }
            $have_contigs = 1;
        }
        close $LINK;

        if($have_contigs) {
            print $OUTPUT2 "$_\n";
            push @contig_2to4, $lines_links[0];
        }else{
            print $OUTPUT1 "$_\n";
        }
    }

    close $INPUT;

    my @contig_3;
    foreach my $contig_2to4_em (@contig_2to4) {
        system "awk '{print \$6}' ./4/$contig_2to4_em.link | sort - | uniq -c - | sed 's/^[ ]*//' >./4/$contig_2to4_em.chr"
            and die "Error: Command failed to run normally: $?\n";

        open my $CHR_FILE, '<', "./4/$contig_2to4_em.chr"
            or die "Error: Cannot create file '$output_dir/link_new/4/$contig_2to4_em.chr': $!\n";
        my $chr_sum = 0;
        my @chr;
        my @chr_num;
        while(<$CHR_FILE>) {
            chomp;
            my @lines_chr = split /\s+/, $_;
            push @chr_num, $lines_chr[0];
            push @chr, $lines_chr[1];
            $chr_sum += $lines_chr[0];
        }
        close $CHR_FILE;

        while (@chr_num) {
            my $chr_num = shift @chr_num;
            my $chr = shift @chr;

            my $pr = $chr_num / $chr_sum;
            if ( $pr >= 0.95 ) {
                open my $CHR_OUT, '>>', "./3.contig.name"
                    or die "Error: Cannot create file '$output_dir/link_new/3.contig.name': $!\n";
                print $CHR_OUT "$contig_2to4_em $chr\n";
                push @contig_3, $contig_2to4_em;
                close $CHR_OUT;

                move "./4/$contig_2to4_em.link", "./3/$contig_2to4_em.link"
                    or die "Error: Cannot move file '$output_dir/link_new/4/$contig_2to4_em.link': $!\n";
                unlink "./4/$contig_2to4_em.chr";
                last;
            }
        }
        unlink "./4/$contig_2to4_em.chr" if ( -e "./4/$contig_2to4_em.chr" );
    }

    mkdir "./2/2a";
    mkpath "./2/2b/left";
    mkpath "./2/2b/right";

    foreach my $contig_3_em ( @contig_3 ) {
        system "sort -u ./3/$contig_3_em.link >./3/$contig_3_em.link.new"
            and die "Error: Command failed to run normally: $?\n";
    }

    system "awk 'NR==FNR{a[\$1]=\$2;next}{print \$1,a[\$1],\$2}' ./contig.name ./3.contig.name >./3.contig.name.length"
            and die "Error: Command failed to run normally: $?\n";


    open my $CONTIG_LENGTH, '<', "./3.contig.name.length"
        or die "Error: Cannot open file '$output_dir/link_new/3.contig.name.length': $!\n";

    while (<$CONTIG_LENGTH>) {
        ( my $contig3_contig, my $contig3_length, my $contig3_chr ) = split /\s+/, $_;
        do3($contig3_contig, "./3/$contig3_contig.link.new", $contig3_length, $contig3_chr);
    }

    close $CONTIG_LENGTH;

    system 'for i in `ls ./4/`; do if [ -s ./4/$i ]; then echo $i | cut -d . -f1 >>./4.contig.name;fi; done'
        and die "Error: Command failed to run normally: $?\n";

    system "awk 'NR==FNR{a[\$1]=\$2;next}{print \$1,a[\$1],\$2}' ./contig.name ./4.contig.name > ./4.contig.name.length"
            and die "Error: Command failed to run normally: $?\n";

    open $CONTIG_LENGTH, '<', "./4.contig.name.length"
        or die "Error: Cannot open file '$output_dir/link_new/4.contig.name.length': $!\n";

    while (<$CONTIG_LENGTH>) {
        ( my $contig4_contig, my $contig4_length ) = split /\s+/, $_;
        do4($contig4_contig, "./4/$contig4_contig.link", $contig4_length);
    }

    close $CONTIG_LENGTH;

    system 'for i in `ls ./4/`; do if [ -s ./4/$i ]; then echo $i | cut -d . -f1 >>./1/4.name;fi; done'
        and die "Error: Command failed to run normally: $?\n";

    system 'for i in `ls ./3/`; do if [ -s ./3/$i ]; then echo $i | cut -d . -f1 >>./1/3.name;fi; done'
        and die "Error: Command failed to run normally: $?\n";

    open my $CONTIG_NAME, '<', "./1/3.name"
        or die "Error: Cannot open file '$output_dir/link_new/1/3.name': $!\n";

    while (<$CONTIG_NAME>) {
        chomp;
        system "awk '{print \$6}' ./3/$_.link | sort - | uniq -c - | sed 's/^[ ]*//' - | sort -nrk 1 - | head -n 1 - > ./3/$_.chr"
            and die "Error: Command failed to run normally: $?\n";
        system "chrpart3=$(awk '{print \$2}' ./3/$_.chr); printf \"$_ \$chrpart3\\n\" >> ./1/3.name.re"
            and die "Error: Command failed to run normally: $?\n";
    }

    unlink "./1/2to4.name";

    system "awk 'NR==FNR{a[\$1]=\$2;next}{print \$1,a[\$1],\$2}' ./contig.name ./1/4.name > ./1/4.name.re"
        and die "Error: Command failed to run normally: $?\n";

    system 'mv ./1/4.name.re ./1/4.name'
        and die "Error: Command failed to run normally: $?\n";
    system 'mv ./1/3.name.re ./1/3.name'
        and die "Error: Command failed to run normally: $?\n";

    system 'for i in `ls ./4/`; do if [ ! -s ./4/$i ]; then echo $i | cut -d . -f1 >>./1/5.name;fi; done'
        and die "Error: Command failed to run normally: $?\n";

    system "rm -rf ./4";
    mkdir "./4";

    return 1;
}

sub do3 {
    ( my $opt_n, my $opt_i, my $opt_l, my $opt_c) = @_;

    open IN,"<$opt_i";
    open OUT1,">$opt_i.left";
    open OUT2,">$opt_i.right";
    while(<IN>) {
        chomp $_;
        my @a = split /\s+/, $_;
        my $ed = $opt_l - 500;
        if( ($a[4] <= 500) && ($a[3] <= 500) && $a[5] eq $opt_c) {
            print OUT1 "$_\n"; ##output the left end results
        } elsif(($a[4] >= $ed) && ($a[3] >= $ed) && $a[5] eq $opt_c) {
            print OUT2 "$_\n"; ##output the right end results
        }
    }
    close IN;
    close OUT1;
    close OUT2;

    my $left = 0;
    my $right = 0;
    my %hash;  # chromosome alignment results

    # left end place
    foreach my $file_suffix ( qw \ left right\ ){
        if ( -s "$opt_i.$file_suffix" ) {
            open IN,"<$opt_i.$file_suffix";
            my @st;
            my @ed;
            while(<IN>) {
                my @cc=split/\s+/,$_;
                ($cc[7] , $cc[8]) = ($cc[8] , $cc[7]) if ( $cc[7] > $cc[8] );
                push @st, $cc[7];
                push @ed, $cc[8];
            }
            close IN;

            if( max( @ed ) - min( @st ) <= 2000 ) {
                $left  = 1 if $file_suffix eq "left";   # left  ensure
                $right = 1 if $file_suffix eq "right";  # right ensure
                my $stmid = int(mid(@st));
                my $edmid = int(mid(@ed));
                $hash{$opt_n}{$file_suffix}="$opt_n $stmid $edmid";  ##left end on chromosome
            }
        }
    }

    open OUT1,">>./1/2bleft.name";
    open OUT2,">>./1/2bright.name";
    open OUT3,">>./1/2a.name";

    if($left==0 && $right==0) {
        unlink "$opt_i","$opt_i.right","$opt_i.left";
    }
    elsif($left != 0 && $right==0) {
        my @a=split/\s+/,$hash{$opt_n}{"left"};
        print OUT1 "$opt_n $opt_l $opt_c $a[1] $a[2]\n";
        system "mv ./3/$opt_n.link ./2/2b/left";
        unlink "$opt_i","$opt_i.right","$opt_i.left";
    }
    elsif($left == 0 && $right != 0) {
        my @a=split/\s+/,$hash{$opt_n}{"right"};
        print OUT2 "$opt_n $opt_l $opt_c $a[1] $a[2]\n";
        system"mv ./3/$opt_n.link ./2/2b/right";
        unlink "$opt_i","$opt_i.right","$opt_i.left";
    }
    elsif($left != 0 && $right != 0) {
        my @a=split/\s+/,$hash{$opt_n}{"left"};
        my @b=split/\s+/,$hash{$opt_n}{"right"};
        print OUT3 "$opt_n $opt_l $opt_c $a[1] $a[2] $b[1] $b[2]\n";
        system"mv ./3/$opt_n.link ./2/2a";
        unlink "$opt_i","$opt_i.right","$opt_i.left";
    }

    sub mid{
        my @list = sort{$a<=>$b} @_;
        my $count = @list;
        return undef if( $count == 0 );

        if( ($count%2) == 1 ){
            return $list[int(($count-1)/2)];
        } elsif ( ($count%2)==0 ) {
            return ($list[int(($count-1)/2)]+$list[int(($count)/2)])/2;
        }
    }

    return 1;
}

sub do4 {
    ( my $opt_n, my $opt_i, my $opt_l) = @_;

    open IN,"<$opt_i";
    open OUT1,">$opt_i.left";
    open OUT2,">$opt_i.right";
    while(<IN>) {
        chomp;
        my @a=split/\s+/,$_;
        my $ed=$opt_l-500;
        if(($a[4] <= 500) && ($a[3] <= 500)) {
            print OUT1 "$_\n"; # output the left end results
        } elsif(($a[4] >= $ed) && ($a[3] >= $ed)) {
            print OUT2 "$_\n"; ##output the right end results
        }
    }
    close IN;
    close OUT1;
    close OUT2;

    ##left end place
    PLACE:foreach my $file_suffix ( qw \ left right\ ){
        system"awk '{print \$6}' $opt_i.$file_suffix | sort - | uniq -c - | sed 's/^[ ]*//' > $opt_i.$file_suffix.chr";
        if ( -s "$opt_i.$file_suffix.chr") {
            open IN,"<$opt_i.$file_suffix.chr";
            my $sum=0;
            my @b;
            my @chr;
            while(<IN>) {
                chomp;
                my @a = split /\s+/, $_;
                push @b, $a[0];
                push @chr, $a[1];
                $sum += $a[0];
            }
            close IN;
            foreach my $b_em ( @b ) {
                if( $b_em/$sum >=0.95) {
                    move "$opt_i", "./3";
                    unlink "$opt_i.$file_suffix", "$opt_i.$file_suffix.chr";
                    last PLACE;
                }
            }
        }
        unlink "$opt_i.$file_suffix", "$opt_i.$file_suffix.chr";
    }
    unlink "$opt_i.right" if ( -e "$opt_i.right" );

    return 1;
}

1;

__END__
