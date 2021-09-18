#!/usr/bin/perl
open IN,"<$ARGV[0]";
open OUT,">./1/1.name";
open OUT1,">./1/2to4.name";
while(<IN>)
{
	chomp $_;
	my @b=split/\s+/,$_;
	my $n=0;
	open IN1,"<$ARGV[1]";
	open OUT2,">./4/$b[0].link";
	while(my $ll=<IN1>)
	{
		my @a=split/\s+/,$ll;
		if($b[0] eq $a[1])
		{
			print OUT2 $ll;
			$n++;
		}
	}
	close IN1;
	if($n==0)
	{
		print OUT "$_\n";
	}
	else
	{
		print OUT1 "$_\n";
	}
}
close IN;