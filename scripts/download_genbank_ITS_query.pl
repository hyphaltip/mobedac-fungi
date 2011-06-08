#!/usr/bin/perl -w
use strict;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Getopt::Long;

my $basedir = 'ITSdb';
my $mindate = '09/01/2010';
my $query_string = '(internal transcribed spacer) AND Fungi[Organism]';
my $gb = Bio::DB::GenBank->new();

my $query;
if( $mindate ) {
   $query = Bio::DB::Query::GenBank->new(-db=>'nucleotide',
                                         -query =>$query_string,
					 -mindate=>$mindate,
                                         );
} else {
   $query = Bio::DB::Query::GenBank->new(-db=>'nucleotide',
                                         -query=>$query_string,
                                         );
}
mkdir($basedir) unless -d $basedir;
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $date = "$year-$mon-$mday";
my $out = Bio::SeqIO->new(-file => ">$basedir/$date.gbk",-format=>'genbank');
my $stream = $gb->get_Stream_by_query($query);
while (my $seq = $stream->next_seq) {
 $out->write_seq($seq);
}
