#!/usr/bin/perl -w
use strict;
use warnings;
use Bio::DB::Taxonomy;
use File::Spec;
use Bio::SeqIO;

my $nodesfile = "./ncbi/nodes.dmp";
my $namefile = "./ncbi/names.dmp";
my $debug = 0;
my $taxdb = Bio::DB::Taxonomy->new(-source => 'flatfile',
                                   -nodesfile => $nodesfile,
                                   -namesfile => $namefile);
warn( "loaded taxonomy db..done\n") if $debug;

my $seqfile = shift;

my $in = Bio::SeqIO->new(-format => 'fasta', -file => $seqfile);
my $out = Bio::SeqIO->new(-format => 'fasta');
while(my $seq = $in->next_seq ) {
    my $desc = $seq->description;
    my (@hierarchy) = split(/;/,$desc);
    my $taxon_name = pop @hierarchy;
    my $id = $taxdb->get_taxonid($taxon_name);
    $seq->description("$id $desc");
    $out->write_seq($seq);
}
