#!/usr/bin/perl
use Getopt::Std;
use List::Util qw/max min/;

use vars qw($opt_n $opt_i $opt_l $opt_c);
getopts('n:i:l:c:');

open IN,"<$opt_i";
open OUT1,">$opt_i.left";
open OUT2,">$opt_i.right";
while(<IN>)
{
	chomp $_;
	my @a=split/\s+/,$_;
	my $ed=$opt_l-500;
	if(($a[4] <= 500) && ($a[3] <= 500) && $a[5] eq $opt_c)
	{
		print OUT1 "$_\n"; ##output the left end results
		$n++;
	}
	elsif(($a[4] >= $ed) && ($a[3] >= $ed) && $a[5] eq $opt_c)
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

##left end place
if (-z "$opt_i.left")
{
	;
}
else
{
	open IN,"<$opt_i.left";
	my @st,@ed;
	my $num=0;
	while(<IN>)
	{
		my @cc=split/\s+/,$_;
		$st[$num]=$cc[7];
		$ed[$num]=$cc[8];
		$num++;
	}
	close IN;
	for $i(0..$num-1)
	{
		if($st[$i] > $ed[$i])
		{
			$it=$st[$i];
			$st[$i]=$ed[$i];
			$ed[$i]=$it;
		}
	}
	my $stmin=min(@st);
	my $edmax=max(@ed);
	my $substr=$edmax-$stmin;
	if($substr <= 2000)
	{
		$left=1;   ##left ensure
		my $stmid=int(mid(@st));
		my $edmid=int(mid(@ed));
		$hash{$opt_n}{"left"}="$opt_n $stmid $edmid";  ##left end on chromosome
	}
	undef @st;
	undef @ed;
}

##right end place
if (-z "$opt_i.right")
{
	;
}
else
{
	open IN,"<$opt_i.right";
	my @st,@ed;
	my $num=0;
	while(<IN>)
	{
		my @cc=split/\s+/,$_;
		$st[$num]=$cc[7];
		$ed[$num]=$cc[8];
		$num++;
	}
	close IN;
	for $i(0..$num-1)
	{
		if($st[$i] > $ed[$i])
		{
			$it=$st[$i];
			$st[$i]=$ed[$i];
			$ed[$i]=$it;
		}
	}
	my $stmin=min(@st);
	my $edmax=max(@ed);
	my $substr=$edmax-$stmin;
	if($substr <= 2000)
	{
		$right=1;  ##right ensure
		my $stmid=int(mid(@st));
		my $edmid=int(mid(@ed));
		$hash{$opt_n}{"right"}="$opt_n $stmid $edmid";  ##right end on chromosome
	}
	undef @st;
	undef @ed;
}


open OUT1,">>./1/2bleft.name";
open OUT2,">>./1/2bright.name";
open OUT3,">>./1/2a.name";
if($left==0 && $right==0)
{
	unlink "$opt_i","$opt_i.right","$opt_i.left";
}
elsif($left ne 0 && $right==0)
{
	my @a=split/\s+/,$hash{$opt_n}{"left"};
	print OUT1 "$opt_n $opt_l $opt_c $a[1] $a[2]\n";
	system"mv ./3/$opt_n.link ./2/2b/left";
	unlink "$opt_i","$opt_i.right","$opt_i.left";
}
elsif($left == 0 && $right ne 0)
{
	my @a=split/\s+/,$hash{$opt_n}{"right"};
	print OUT2 "$opt_n $opt_l $opt_c $a[1] $a[2]\n";
	system"mv ./3/$opt_n.link ./2/2b/right";
	unlink "$opt_i","$opt_i.right","$opt_i.left";
}
elsif($left ne 0 && $right ne 0)
{
	my @a=split/\s+/,$hash{$opt_n}{"left"};
	my @b=split/\s+/,$hash{$opt_n}{"right"};
	print OUT3 "$opt_n $opt_l $opt_c $a[1] $a[2] $b[1] $b[2]\n";
	system"mv ./3/$opt_n.link ./2/2a";
	unlink "$opt_i","$opt_i.right","$opt_i.left";
}

sub mid{
    my @list = sort{$a<=>$b} @_;
    my $count = @list;
    if( $count == 0 )
    {
        return undef;
    }
    if(($count%2)==1){
        return $list[int(($count-1)/2)];
    }
    elsif(($count%2)==0){
        return ($list[int(($count-1)/2)]+$list[int(($count)/2)])/2;
    }
}
