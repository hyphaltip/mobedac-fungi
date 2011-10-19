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
use TokyoCabinet; # this is much faster than DB_File (BerkeleyDB) for millions of key-value pairs
use Time::HiRes qw(gettimeofday tv_interval);

my $acc_count = 247143112;
my $gi_count  = 230829154;
my $chunk = 1_000_000;
my $gbacc2gi = '/scratch/gbacc/GbAccList.gz';
my $gi2taxon = '/scratch/gbacc/gi_taxid_nucl.dmp.gz';
my $acc2id_idx = '/scratch/gbacc/acc2gi.tch';
my $gi2taxa_idx = '/scratch/gbacc/gi2taxon.tch';
my $nodesfile = "./ncbi/nodes.dmp";
my $namefile = "./ncbi/names.dmp";

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


## Make local taxonomy db
warn( "Constructing flatfile db..\n") if $debug;
my $taxdb = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                   -nodesfile => $nodesfile,
                                   -namesfile => $namefile);
warn( "loaded taxonomy db..done\n") if $debug;

# create the object
my $acc2gi_hdb = TokyoCabinet::HDB->new();


# open the ACC2GI database    
my $ready = 0;
if( -f $acc2id_idx && ! -z $acc2id_idx ) {
    # this code opens the DB 
    if(!$acc2gi_hdb->open($acc2id_idx)){
	my $ecode = $acc2gi_hdb->ecode();
	die ("open error: %s\n", $acc2gi_hdb->errmsg($ecode));
    }
    
    if( $acc2gi_hdb->rnum == 0 ) {
	# if it exists but is empty, close and reopen
	$acc2gi_hdb->close;
	unlink($acc2id_idx);
	$ready = 0;
    } else { $ready = 1; }
}
unless( $ready ) {
    $acc2gi_hdb->tune(int($acc_count*2),undef,undef,$acc2gi_hdb->TLARGE | 
		      $acc2gi_hdb->TTCBS); # 1 * num of Accs
    if(!$acc2gi_hdb->open($acc2id_idx, $acc2gi_hdb->OWRITER | $acc2gi_hdb->OCREAT)){
	my $ecode = $acc2gi_hdb->ecode();
	die ("open error: %s\n", $acc2gi_hdb->errmsg($ecode));
    }
    open(my $accfh => "zcat $gbacc2gi |") || die $!;
    if( $acc2gi_hdb->rnum == 0 ) {
	my $i = 0;
	my $lasttime = [gettimeofday];
	while(<$accfh>) {
	    chomp;
	    my ($acc,$ver,$gi) = split(/,/,$_);
# store records
	    if(!$acc2gi_hdb->put($acc, $gi) ){
		my $ecode = $acc2gi_hdb->ecode();
		printf STDERR ("put error: %s\n", $acc2gi_hdb->errmsg($ecode));
		
	    }
	    if( $i++ % $chunk == 0 && $i > 1) { 
		warn(sprintf("$i %d seconds for $chunk records\n", tv_interval($lasttime))) if $debug;
		$lasttime = [gettimeofday];
	    }
	}
	$acc2gi_hdb->sync();
	#$acc2gi_hdb->optimize();
    }
}
warn($acc2gi_hdb->rnum, " acc2gi elements\n") if $debug;

# open the GI2TAXA database    
my $gi2taxa_hdb = TokyoCabinet::HDB->new();
$ready = 0;
if( -f $gi2taxa_idx && ! -z $gi2taxa_idx ) {
    if(!$gi2taxa_hdb->open($gi2taxa_idx)){
	my $ecode = $gi2taxa_hdb->ecode();
	die ("open error: %s\n", $gi2taxa_hdb->errmsg($ecode));
    }
    if( $gi2taxa_hdb->rnum == 0 ) {
	# if it exists but is empty, close and reopen
	$gi2taxa_hdb->close;
	unlink($gi2taxa_idx);
	$ready = 0;
    } else {
	$ready = 1;
    }   
} 

unless(  $ready ) {
    $gi2taxa_hdb->tune(int($gi_count*2),undef,undef,
		       $gi2taxa_hdb->TLARGE | $gi2taxa_hdb->TTCBS); # 1 * num of gi2taxa
    # open the database    
    if(! $gi2taxa_hdb->open($gi2taxa_idx, $gi2taxa_hdb->OWRITER | $gi2taxa_hdb->OCREAT)){
	my $ecode = $gi2taxa_hdb->ecode();
	die ("open error: %s\n", $gi2taxa_hdb->errmsg($ecode));
    }
    
    open(my $gitaxafh => "zcat $gi2taxon |") || die $!;    
    if( $gi2taxa_hdb->rnum == 0 ) {
	my $i = 0;
	my $lasttime = [gettimeofday];
	while(<$gitaxafh>) {
	    my ($gi,$taxid) = split;
	    # store records
	    if(!$gi2taxa_hdb->put($gi, $taxid) ){
		my $ecode = $gi2taxa_hdb->ecode();
		printf STDERR ("put error: %s\n", $gi2taxa_hdb->errmsg($ecode));
	    }
	    if( $i++ % $chunk == 0 && $i > 1) { 
		warn(sprintf("$i %d seconds $chunk records\n",tv_interval($lasttime))) if $debug;
		$lasttime = [gettimeofday];
	    }
	}
	$gi2taxa_hdb->sync();
	#$gi2taxa_hdb->optimize();
    }
}
warn($gi2taxa_hdb->rnum, " gi2taxid elements\n") if $debug;


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
while (my $seq = $seqio->next_seq) {
    
    my $line = $seq->display_name;
    #print $line,"\n";
    ## GI/Accession number is the first item before the first underscore
    ## Organism name/Taxonomy ID is everything after the first underscore
    my ($ID,@descriptors) = split(/\_/, $line);
    
    my $gi = $acc2gi_hdb->get($ID);

    my ($input,$taxid,$curr,@hierarchy);
    my $input = join(" ",@descriptors);	

    if( defined $gi && defined($taxid = $gi2taxa_hdb->get($gi))) {
	# found taxid	
	
    } else {
	## If the input provided is a Taxonomy ID number, 
	##  get the corresponding organism name
	if( $input =~ /^\d+$/) {
	    $taxid = $input;
	}	
    }
    if( $taxid ) {
	$curr = $taxdb->get_taxon(-taxonid => $taxid);
	if( ! defined $curr ) { 
	    warn("cannot lookup $taxid should I check for $input instead?\n"); 
	    $curr = $taxdb->get_taxon(-name => $input);
        }
    } else {
	$curr = $taxdb->get_taxon(-name => $input);
    }
    unless( $curr ) { 
	## Flag missing taxonomic ranks
	$unidentified = 1;
	print UNID "$ID $input\n";
	next;
    }
    while( defined $curr ) {
	if( $curr->rank ne 'no rank' || $show ) {
	    unshift(@hierarchy,getName($curr));
	}
	$curr = $curr->ancestor; # go up a level in hierarchy    
    }
    $seq->description(join(";",@hierarchy).".");
    $out->write_seq($seq);      
}
warn( "..done\n") if $debug;

if($unidentified) {
  print "Some organisms were not found in the taxonomy database.\n";
  print "These organisms are in the \"$filename\.unidentified\" file.\n";
}

# close the database
if(!$acc2gi_hdb->close()){
    my $ecode = $acc2gi_hdb->ecode();
    die ("close error: %s\n", $acc2gi_hdb->errmsg($ecode));
}


# close the database
if(!$gi2taxa_hdb->close()){
    my $ecode = $gi2taxa_hdb->ecode();
    die ("close error: %s\n", $gi2taxa_hdb->errmsg($ecode));
}
 


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
