#!/usr/bin/perl -w
# Script: seqcount.pl
# Description: Displays the number of sequences in all fasta (default) or genbank files found in the current directory, or a specific fasta/genbank file
# Author: Steven Ahrendt
# email: sahrendt0@gmail.com
# Date: 9.1.11
#         v1.0
#         v1.5  : added genbank support
#         v1.6  : added file input option
#################################
# Usage: seqcount.pl [type] [output] [file]
#################################
# "type":   -g for genbank
# "output": -v for verbose (create .seqlist/count files)
# "file":   -f, followed by the filename
#################################
# .seqlist files contain the display ids for each sequence, 
#  "counts" file containing the name of each sequence file and the number
#  of sequences within it
#################################

use strict;
use Bio::Seq;
use Bio::SeqIO;

my $format = "fasta";
my $ext = "f.*a";
my $verbose = 0;
my @files;

# Handle command-line options
use Getopt::Std;
my %opts;
getopts('gvf', \%opts);
if (exists $opts{'g'}) 
{
  $format = "genbank";
  $ext = "gbk";
}
if (exists $opts{'v'}) 
{
  $verbose = 1;
}
if (exists $opts{'f'}) 
{
  push(@files,$ARGV[0]);
}
else
{
  my $dir = ".";
  opendir(DIR,$dir);
  @files = sort(grep {/\.$ext$/} readdir(DIR));
  closedir(DIR);
}
my $numfiles = @files;

if($numfiles > 0)
{
  if($verbose){open(COUNTS,">counts");}
  foreach my $file (@files)
  {
    my $numseqs=0;
    my $seqfile = Bio::SeqIO->new(-file=>$file,
                                    -format=>$format);
    if($verbose){open(LIST,">$file.seqlist");}
    while(my $seq = $seqfile->next_seq())
    {
	    $numseqs++;
	    if($verbose){print LIST $seq->display_id(),"\n";}
    }
    if($verbose){close(LIST);}
    print "$file\t$numseqs\n";
    if($verbose){print COUNTS "$file\t$numseqs\n";}
  }
  if($verbose){close(COUNTS);}
}
else
{
  print "Error: no $format files found in current directory.\n";
}
