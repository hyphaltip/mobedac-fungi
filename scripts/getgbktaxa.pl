#!/usr/bin/perl -w
# Script: getgbktaxa.pl
# Description: Get all unique taxIDs in gbk files
# Author: Steven Ahrendt
# email: sahrendt0@gmail.com
#        sahre001@ucr.edu
# Date: 8.31.11
#          v.1.0  : 
#	   v.1.1  : Updated taxID extraction (sometimes "db_xref" fields are not "taxon" fields)
#          v.1.2  : produces .taxlist file for use with taxonomy.pl
#####################################
# Usage: getgbktaxa.pl [gbkfile]
#####################################
# If a gbkfile is not provided, it will gather all .gbk files in ITSdb directory
#####################################
use strict;
use Bio::SeqIO;

my $dir = "ITSdb";

my @gbkfiles;
if ($ARGV[0]) 
{
  push(@gbkfiles, $ARGV[0]);
}
else
{
  opendir(DIR,$dir);
  @gbkfiles = grep { /\.gbk$/ } readdir(DIR);
  close(DIR);
  @gbkfiles = sort @gbkfiles;
}


foreach my $gbkfile (@gbkfiles)
{
  #print "$gbkfile:\n";
  my $gbk = Bio::SeqIO->new(-file   => "$dir/$gbkfile", 
                            -format => 'genbank');
  my @id_ar;
  open(TAXLIST,">$gbkfile\.taxlist");
  while(my $seq = $gbk->next_seq)
  {
    my @taxonomy = taxonomy($seq);
    my $taxID = $taxonomy[0];
    print TAXLIST $seq->accession();
    print TAXLIST "\_$taxID\n";
    push(@id_ar,$taxID);
  }
  close(TAXLIST);
  my %id_hash = map{$_,1} @id_ar;
  open(SEQLIST, ">$gbkfile\.seqlist");
  foreach my $id (keys %id_hash)
  {
    print SEQLIST $id,"\n";
  }
  close(SEQLIST);
}
#print `perl getchytrid.pl`;


## Subroutine to extract taxID and species name
#  Author: Steven Ahrendt
sub taxonomy {
   my ($obj,@args) = @_;
   my @t;
   my @source_features = grep {$_->primary_tag eq 'source' } $obj->get_SeqFeatures;
   my @dbxref = (($source_features[0])->get_tag_values('db_xref'));
   foreach my $dbxref (@dbxref)
   {
     if($dbxref =~ m/taxon/)
     {
       push @t, (split(":",$dbxref))[1]; # "taxon:taxID"
       push @t, $source_features[0]->get_tag_values('organism');
     }
   }
   return @t;
} # taxonomy() subroutine
