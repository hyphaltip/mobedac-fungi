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
use Bio::DB::Taxonomy;
use File::Spec;
use Bio::SeqIO;
use Getopt::Long;

my ($name);
my $numseqs = 0;
my $unidentified = 0; # Unidentified organisms
my $show = 0; # Show or hide ranks with value "no rank"
my $outext = "";
# Handle command-line options
my $debug;

GetOptions( 's|show!'=> sub { $show = 1; $outext = '.norank'},
	    'verbose|debug!' => \$debug,
	    );

my $IDfile = shift @ARGV;


## Function to process a scientific name based on taxonomic rank
#   Returns: the name in angle brackets "<>" if the rank is "no rank"
sub getName {
  my $to = shift;
  my $rank = $to->rank();
  my $res = $to->scientific_name;
  if($rank eq "no rank") {
      $res = join("","<",$to->scientific_name,">");
  }
  return $res;
}

## Make local taxonomy db
warn( "Constructing flatfile db..\n") if $debug;
my $nodesfile = "./ncbi/nodes.dmp";
my $namefile = "./ncbi/names.dmp";
my $taxdb = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                   -nodesfile => $nodesfile,
                                   -namesfile => $namefile);
warn( "loaded taxonomy db..done\n") if $debug;

## Store filename for output files
my @f = split(/\./,$IDfile);
pop(@f);
my $filename = join(".",@f);
$filename = join("",$filename,$outext);

## Get descriptions from ID file
my $seqio = Bio::SeqIO->new(-format => 'fasta',
			    -file   =>  $IDfile);
warn "Creating taxonomy hierarchy..\n" if $debug;

my $out = Bio::SeqIO->new(-format=> 'fasta',
			  -file  => ">$filename\.taxonomy.seqs");
open(UNID,">$filename\.unidentified");
open(OUT,">$filename\.taxonomy");
while (my $seq = $seqio->next_seq) { 
    my $line = $seq->display_name;
    #print $line,"\n";
    ## GI/Accession number is the first item before the first underscore
    ## Organism name/Taxonomy ID is everything after the first underscore
    my ($ID,@descriptors) = split(/\_/, $line);
    my $org = join(" ",@descriptors);
    
    ## If the input provided is a Taxonomy ID number, 
    ##  get the corresponding organism name
    my $input = $org;
    if($input =~ /^\d+$/) {
	$input = ($taxdb->get_taxon(-taxonid => $org))->scientific_name;
    }

  ## Using organism name, generate hierarchy
  my @hierarchy;
  my %ranks;
  if(my $curr = $taxdb->get_taxon(-name => $input)) {
      my $name = getName($curr);
      #print $name,";";
      #unshift(@hierarchy,$name);
      ## Default: hide "no rank"
      ## Flag: -s = show "no rank"

      # could also just check and see if $curr->rank is NULL?
      while($curr) {
	  
          if($curr->rank eq "kingdom" || $curr->rank eq "phylum" || $curr->rank eq "class" || $curr->rank eq "order" || $curr->rank eq "family" || $curr->rank eq "genus")
          {
	  #print "<",$curr->rank,">;";
	  $name = getName($curr);
          #print $name,";";
	  if($name !~ /^</ || $show)  {
	      unshift(@hierarchy,$name);
	  }
          }
	  $curr = $curr->ancestor;
      }
      #print "\n";
      print OUT $ID,"\t",join(";",@hierarchy),"\n";
      $seq->description(join(";",@hierarchy).".");
      $out->write_seq($seq);      
  } else {
      ## Flag missing taxonomic ranks
      $unidentified = 1;
      print UNID "$ID $input\n";
  }

}
warn( "..done\n") if $debug;

if($unidentified)
{
  print "Some organisms were not found in the taxonomy database.\n";
  print "These organisms are in the \"$filename\.unidentified\" file.\n";
}
close(UNID);
close(OUT);
