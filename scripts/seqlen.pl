#!/usr/bin/perl -w
# Script: seqlen.pl
# Description: Extracts sequence lengths from a fasta file
# Author: Steven Ahrendt
# email: sahrendt0@gmail.com
# Date: 6.17.11
#	v 1.0:	Can probably make more general/robust
#		Make a simple histogram
#		Open data in oocalc
#	v 1.5:	Updated w/ bioperl
#		No histogram
#		No oocalc
#		works for any fastafile
#	v 1.6:	added stats
#	v 1.7:  correct "div-by-0" error
#####################################################
# Usage: seqlen.pl fastafile
#####################################################

use strict;
use Bio::Perl;
use Bio::SeqIO;
use List::Util qw[min max];

sub stddev {
	my @a = @_;
	my $n = @a;
	my $c = (sum(@a) ** 2)/$n;
	my $num = sumsqr(@a) - $c;

	#print "\n";
	#print sumsqr(@a);
	#print " - ";
	#print "(",sum(@a),"^2)/($n)\n";
	#print "-----------------------\n";
	#print $n,"-1\n";

	return ($num/($n-1)) ** 0.5;
}

sub mean {
	my @a = @_;
	my $n = @a;
	return (sum(@a)/$n);
}


sub sumsqr {
	my @a = @_;
	my $ret = 0;
	foreach my $i (@a)
	{
		$ret += ($i**2);
	}
	return $ret;
}

sub sum {
	my @a = @_;
	my $ret = 0;
	foreach my $i (@a)
	{
		$ret += $i;
	}
	return $ret;
}

sub round {
	my($number) = shift;
	return int($number + .5);
}

sub infilename {
	my $file = shift;
	return substr($file,0,index($file,"."));
}

sub histogram{
	my @data = @_;
	my $size = @data;
	my $bins = 10;
	my $min = 0;
	my $max = max(@data);
	print "$max\n";
	if(($max%10) != 0)
	{
		$max /= 10;
		$max = round($max);
		$max *= 10;
	}
	my $c = $min;
	print "|";
	while($c<$max)
	{
		printf(" %3d |",$c);
		$c += ($max/10);
	}
	print " $max |";
	print "\n";
}

my (@tmp,@len);
my @lengths;
my $infile = $ARGV[0];
print `perl -pi -e 's/\r/\n/g' $infile`;

my $numseqs = 0;

my $fastafile = Bio::SeqIO->new(-file=>"<$infile", -format=>"fasta");
my $csvfile = join(".",infilename($infile),"csv");

open(OUT,">$csvfile");
while(my $seq = $fastafile->next_seq)
{
	push(@lengths,$seq->length);
	print OUT $seq->display_id,",",$seq->length,"\n";
}
close(OUT);

$numseqs = @lengths;
print "# Seqs:\t$numseqs\n";
print "Min:\t",min(@lengths),"\n";
print "Max:\t",max(@lengths),"\n";
printf("Mean:\t%.2f\n",mean(@lengths));
if($numseqs == 1)
{
	printf("Stdv:\tNA\n");
}
else
{
	printf("Stdv:\t%.2f\n",stddev(@lengths));
}
print "-----\n";
print "Lengths are in \"$csvfile\"\n";
print "-----\n";
