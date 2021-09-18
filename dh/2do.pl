#!/usr/bin/perl
open IN,"<$ARGV[0]";
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
		open OUT,">>./3.contig.name";
		print OUT "$ARGV[1] $chr[$i]\n";
		close OUT;
		system"mv ./4/$ARGV[1].link ./3";
		unlink"./4/$ARGV[1].chr";
		last;
	}
}
unlink"./4/$ARGV[1].chr";
undef @b,@chr;
