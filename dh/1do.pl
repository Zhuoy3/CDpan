#!/usr/bin/perl
my %hash;
my $line=0;
open IN1,"<$ARGV[1]";
while(my $ll=<IN1>)
{
	chomp $ll;
	my @a=split/\s+/,$ll;
	$hash{$a[1]}{$line}=$ll;
	$line++;
}
close IN1;
open IN,"<$ARGV[0]";
open OUT,">./1/1.name";
open OUT1,">./1/2to4.name";
while(<IN>)
{
	chomp $_;
	my @b=split/\s+/,$_;
	my $n=0;
	open OUT2,">./4/$b[0].link";
	if(exists $hash{$b[0]})
	{
		foreach my $k2 (keys %{$hash{$b[0]}})
		{
			print OUT2 "$hash{$b[0]}{$k2}\n";
		}
		$n++;
	}
	close OUT2;
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
