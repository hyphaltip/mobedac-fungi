#!/usr/bin/perl -w
# Script: seqcount.pl
# Description: Displays the number of sequences in the fasta (default) or genbank files found in the current directory
# Author: Steven Ahrendt
# email: sahrendt0@gmail.com
# Date: 6.16.11
#         v1.0
#         v1.5  : added genbank support
#################################
# Usage: seqcount.pl [type]
#################################
# "type": -g for genbank
#################################
# The script also creates a *.seqlist file containing the ids for each sequence, 
#  as well as a "counts" file containing the name of each sequence file and the number
#  of sequences within it
#################################

use strict;
use Bio::Seq;
use Bio::SeqIO;

my $format = "fasta";
my $ext = "f";

# Handle command-line options
use Getopt::Std;
my %opts;
getopts('g', \%opts);
if (exists $opts{'g'}) 
{
  $format = "genbank";
  $ext = "g";
}

my $dir = ".";

opendir(DIR,$dir);
my @files;
@files = grep {/\.$ext.+$/} readdir(DIR);
closedir(DIR);

my $numfiles = @files;

if($numfiles > 0)
{
  open(COUNTS,">counts");
  foreach my $file (@files)
  {
    my $numseqs=0;
    my $seqfile = Bio::SeqIO->new(-file=>$file,
                                    -format=>$format);
    open(LIST,">$file.seqlist");
    while(my $seq = $seqfile->next_seq())
    {
	    $numseqs++;
	    print LIST $seq->display_id(),"\n";
    }
    close(LIST);
    print "$file\t$numseqs\n";
    print COUNTS "$file\t$numseqs\n";
  }
  close(COUNTS);
}
else
{
  print "Error: no $format files found in current directory.\n";
}
