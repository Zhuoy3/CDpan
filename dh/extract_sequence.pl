#!/usr/bin/perl
my %id;
open IN,"<$ARGV[0]";
while(<IN>)
{
	chomp $_;
	$id{$_}=$_;
}
close IN;

open IN1,"<$ARGV[1]";
open OUT,">./filtered.fa";
my $n=0;
while(my $ll=<IN1>)
{
	chomp $ll;
	if($ll=~/^>/)
	{
		if(exists($id{$ll}))
		{
			print OUT "$ll\n";
			$n=1;
		}
		else
		{
			$n=0;
		}
	}
	else
	{
		if($n==1)
		{
			print OUT "$ll\n";
		}
		else
		{;}
	}
}
close IN1;
close OUT;