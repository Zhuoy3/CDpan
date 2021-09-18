#!/usr/bin/perl
open IN1,"<$ARGV[1]";
my %id;
while(my $ll=<IN1>)
{
	chomp $ll;
	$id{$ll}=$ll;
}
close IN1;

open IN,"<$ARGV[0]";
open OUT,">./keep_id.list";
while(<IN>)
{
	if(/^readID/)
	{;}
	else
	{
		chomp $_;
		my @a=split/\s+/,$_;
		if(exists($id{$a[2]}))
		{
			print OUT ">$a[0]\n";
		}
	}
}
