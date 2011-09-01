#!/usr/bin/perl -w
# Script: taxonomy.pl
# Description: Takes an ID file and determines the taxonomy hierarchy for each species 
# Author: Steven Ahrendt
# email: sahrendt0@gmail.com
#        sahre001@ucr.edu
# Date: 8.30.11
#       v.1.0  :
#               [x] Make local db
#               [x] generate hierarchy
#               [x] identify non-ranked levels
#               [ ] flag missing levels
#####################################
# Usage: taxonomy.pl [-s] IDfile
#####################################
# -s: flag to show or hide ranks with value "no rank"
#     (include flag to show ranks)
#####################################
# ID file should have the format:
#  someID_someTaxID
# where
#  someID = either Accession or GI
# and
#  someTaxID = either NCBI Taxonomy ID or "Genus species" identifier
#####################################

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Bio::DB::Taxonomy;

my ($ID,$org,@lines,$name);
my $numseqs = 0;
my $unidentified = 0; # Unidentified organisms
my $show = 0; # Show or hide ranks with value "no rank"
my $outext = "";
# Handle command-line options
use Getopt::Std;
my %opts;
getopts('s', \%opts);
if (exists $opts{'s'}) 
{
  $show = 1;
  $outext = ".norank";
}

my $IDfile = $ARGV[0];


## Function to process a scientific name based on taxonomic rank
#   Returns: the name in angle brackets "<>" if the rank is "no rank"
sub getName {
  my $to = shift;
  my $rank = $to->rank();
  my $res = $to->scientific_name;
  if($rank eq "no rank")
  {
    $res = join("","<",$to->scientific_name,">");
  }
  return $res;
}

## Make local taxonomy db
print "Constructing flatfile db..\n";
my $nodesfile = "./ncbi/nodes.dmp";
my $namefile = "./ncbi/names.dmp";
my $taxdb = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                   -nodesfile => $nodesfile,
                                   -namesfile => $namefile);
print "..done\n";

## Get descriptions from ID file
open(ID,"<$IDfile") or die "Can't open file \"$IDfile\"\n";
@lines = <ID>;
close(ID);
chomp(@lines);

## Store filename for output files
my @f = split(/\./,$IDfile);
pop(@f);
my $filename = join(".",@f);
$filename = join("",$filename,$outext);

print "Creating taxonomy hierarchy..\n";
open(HASH,">$filename\.taxonomy");
foreach my $line (@lines)
{
  #print $line,"\n";
  ## GI/Accession number is the first item before the first underscore
  ## Organism name/Taxonomy ID is everything after the first underscore
  my @descriptors = split(/\_/, $line);
  $ID = $descriptors[0];
  shift(@descriptors);
  $org = join(" ",@descriptors);

  ## If the input provided is a Taxonomy ID number, 
  ##  get the corresponding organism name
  my $input = $org;
  if($input =~ /^[0-9]+$/)
  {
    $input = ($taxdb->get_taxon(-taxonid => $org))->scientific_name;
  }

  ## Using organism name, generate hierarchy
  my @hierarchy;
  if(my $curr = $taxdb->get_taxon(-name => $input))
  {
    $name = getName($curr);
    #unshift(@hierarchy,$name);
    ## Default: hide "no rank"
    ## Flag: -s = show "no rank"

    if($name !~ /^</)
    {
      unshift(@hierarchy,$name);
    }
    else
    {
      if($show)
      {
        unshift(@hierarchy,$name);
      }
    }
    
    while($curr = $curr->ancestor)
    {
      $name = getName($curr);
      if($name !~ /^</)
      {
        unshift(@hierarchy,$name);
      }
      else
      {
        if($show)
        {
          unshift(@hierarchy,$name);
        }
      }
    }
    print HASH "$ID\t";
    print HASH join(";",@hierarchy);
    print HASH ".\n";
  }
  else
  {
    open(UNID,">>$filename\.unidentified");
    $unidentified = 1;
    print UNID $input,"\n";
    close(UNID);
  }

  ## Flag missing taxonomic ranks
}
print "..done\n";
if($unidentified)
{
  print "Some organisms were not found in the taxonomy database.\n";
  print "These organisms are in the \"$filename\.unidentified\" file.\n";
}
