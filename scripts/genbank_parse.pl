#!/usr/bin/perl -w
# Script: genbank_parser.pl
# Description: Parses a genbank file
# Author: Steven Ahrendt
# email: sahrendt0@gmail.com
#        sahre001@ucr.edu
# Date: 6.17.11
#         v0.5  : lists all taxids and RNA features in a genbank file
#         v1.0  : creates a fasta file with the ITS sequences that are annotated
########################################
# Usage: genbank_parser.pl genbankfile #
########################################

use strict;
use Bio::Seq;
use Bio::SeqIO;

## Subroutine to extract taxID and species name
#  Author: Steven Ahrendt
sub taxonomy {
   my ($obj,@args) = @_;
   my @t;
   my @source_features = grep {$_->primary_tag eq 'source' } $obj->get_SeqFeatures;
   push @t, (split(":",(($source_features[0])->get_tag_values('db_xref'))[0]))[1];
   push @t, $source_features[0]->get_tag_values('organism');
   return @t;
} # taxonomy() subroutine

my $basedir = "ITSdb";
opendir(DIR,"$basedir");
my @gbkfiles = grep { /\.gbk$/ } readdir(DIR);
closedir(DIR);

my $gbkfile = "$basedir/$ARGV[0]";

my $gbk = Bio::SeqIO->new(-file=>$gbkfile, 
                          -format=>'genbank');

my $ITS = Bio::SeqIO->new(-file=>">$ARGV[0]\_ITS.fa",
			       -format=>'fasta');

## Iterate through each sequence in the Genbank file
while(my $seq = $gbk->next_seq)
{
  my $accID = $seq->accession();
  my @taxonomy = taxonomy($seq);
  my $taxID = $taxonomy[0];
  my $species = $taxonomy[1];
  my $total_length = $seq->length();
  print "$taxID\t$accID\t$species\n";
  
  ## Some rRNA features are annotated as 'rRNA'. However some ITS regions
  #   are annotated as 'misc_feature'. The following few lines create one 
  #   array to store all rRNA features
  my @rna_features = grep { $_->primary_tag eq 'rRNA' } $seq->get_SeqFeatures;
  my @misc_features = grep { $_->primary_tag eq 'misc_feature' } $seq->get_SeqFeatures;
  foreach my $mfeat (@misc_features)
  {
    push @rna_features,$mfeat;
  }
  
  ## Iterate through the rRNA features 
  foreach my $feat (@rna_features)
  {
    my @tags = $feat->get_all_tags();
    for my $tag (@tags)
    {
      ## Relevant tags are either 'note', 'product', or 'standard_name'.
      #   As of now, all are taken. Some may be duplicates. (an 'rRNA'
      #   feature can have both a 'note' and a 'product' annotation
      print "     ",$total_length,"\t",$feat->seq->length,"\t",$feat->get_tag_values($tag),"<$tag>","<<",$feat->primary_tag,">>","\n" if (($tag eq 'product') || ($tag eq 'note') || ($tag eq 'standard_name'));
      if(($tag eq 'product') || ($tag eq 'note') || ($tag eq 'standard_name'))
      {
	for my $value ($feat->get_tag_values($tag))
	{	
          if($value =~ m/^I.*T.*S.*$/i)
          {
            ## Write the ITS sequence in fasta format.
            #   Header line will be "taxID|species|description|accession"
            my $ITSobj = Bio::Seq->new(-seq=>$feat->seq->seq,
                                       -desc=>"|$species|$value|$accID",
                                       -id=>$taxID);
            $ITS->write_seq($ITSobj);
          }
        }
      }
    }
  } # For each rna_feature
} # For each genbank sequence



#############################
### UNUSED CODE FRAGMENTS ###
#############################
#  my (%result,@rvalues,$taxID);
#  my $accID = $seq->accession();
#  my $seqlength = $seq->length();
#  push(@rvalues,$accID);
#  for my $feat_object ($seq->get_SeqFeatures) 
#  {  
	#print ":",$accID,"primary tag: ", $feat_object->primary_tag, "\n";
#    if($feat_object->primary_tag() eq "source")
#    {
#      if($feat_object->has_tag("db_xref"))
#      {
#	$taxID = (split(":",($feat_object->get_tag_values("db_xref"))[0]))[1];
        #print $tmp;
        #$taxID = (split(":",$taxID))[1];
	#push @rvalues,$taxID;
#      }
#      push @rvalues, $feat_object->get_tag_values("organism") if ($feat_object->has_tag("organism"));
#    }
#    if(($feat_object->primary_tag() eq "misc_feature") || ($feat_object->primary_tag() eq "rRNA"))
#    {
#      my @tags = $feat_object->get_all_tags();
#      for my $tag (@tags)
#      {
#        push @rvalues, $feat_object->get_tag_values($tag) if (($tag =~ m/product/));# || ($tag =~ m/note/));
	#push @rvalues, "<$tag>";
#      }
#    }
    #if($feat_object->primary_tag() eq "rRNA")
    #{
    #  push @rvalues,$feat_object->get_tag_values("product") if ($feat_object->has_tag("product"));
    #}
    #print $accID,"\t",$taxID,"\t",$species,"\n";
#  }

  #print @rvalues;
#  push @rvalues, $seqlength;

#  $result{$taxID} = \@rvalues;
  #print $result{$accID};
#  for my $key (keys %result)
#  {
#    print "$key:", join(":",@{ $result{$key} }),"\n";
#  }
