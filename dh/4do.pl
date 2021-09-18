#!/usr/bin/perl
use Getopt::Std;
use List::Util qw/max min/;

use vars qw($opt_n $opt_i $opt_l);
getopts('n:i:l:');

open IN,"<$opt_i";
open OUT1,">$opt_i.left";
open OUT2,">$opt_i.right";
while(<IN>)
{
	chomp $_;
	my @a=split/\s+/,$_;
	my $ed=$opt_l-500;
	if(($a[4] <= 500) && ($a[3] <= 500))
	{
		print OUT1 "$_\n"; ##output the left end results
		$n++;
	}
	elsif(($a[4] >= $ed) && ($a[3] >= $ed))
	{
		print OUT2 "$_\n"; ##output the right end results
		$n++;
	}
	else
	{;}
}
close IN;
close OUT1;
close OUT2;

my $left=0;
my $right=0;
my %hash;  ##chromosome alignment results
system"awk '{print $6}' $opt_i.left | sort - | uniq -c - > $opt_i.left.chr";
system"sed -i 's/^[ ]*//' $opt_i.left.chr";
system"awk '{print $6}' $opt_i.right | sort - | uniq -c - > $opt_i.right.chr";
system"sed -i 's/^[ ]*//' $opt_i.right.chr";


##left end place
if (-z "$opt_i.left.chr")
{
	;
}
else
{
	open IN,"<$opt_i.left.chr";
	my $sum=0;
	my $n=0;
	my @b;
	my @chr;
	while(<IN>)
	{
		my @a=split/\s+/,$_;
		$b[$n]=$a[0];
		$chr[$n]=$a[1];
		$sum+=$b[$n];
		$n++;
	}
	close IN;
	for $i(0..$n-1)
	{
		my $pr=$b[$i]/$sum;
		if($pr>=0.95)
		{
			system"mv $opt_i ./3";
			unlink"$opt_i.left","$opt_i.right","$opt_i.left.chr","$opt_i.right.chr";
			goto END;
		}
	}
	undef @b,@chr;
}

##right end place
if (-z "$opt_i.right.chr")
{
	;
}
else
{
	open IN,"<$opt_i.right.chr";
	my $sum=0;
	my $n=0;
	my @b;
	my @chr;
	while(<IN>)
	{
		my @a=split/\s+/,$_;
		$b[$n]=$a[0];
		$chr[$n]=$a[1];
		$sum+=$b[$n];
		$n++;
	}
	close IN;
	for $i(0..$n-1)
	{
		my $pr=$b[$i]/$sum;
		if($pr>=0.95)
		{
			system"mv $opt_i ./3";
			unlink"$opt_i.left","$opt_i.right","$opt_i.left.chr","$opt_i.right.chr";
			last;
		}
	}
	undef @b,@chr;
}

unlink"$opt_i.left","$opt_i.right","$opt_i.left.chr","$opt_i.right.chr";

END:undef @b,@chr;
