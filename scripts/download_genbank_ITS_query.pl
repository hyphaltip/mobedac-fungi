#!/usr/bin/perl -w
use strict;
use Bio::DB::GenBank;
use Bio::DB::Query::GenBank;
use Getopt::Long;

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

my $basedir = 'ITSdb';
my $DEBUG = 0;
my @mindate = (1993,01,01);
my @maxdate = (1900+$year,$mon+1,$mday);
GetOptions(
	'm|min:s' => \$mindate[0],
	'b|base:s' => \$basedir,
	);
my $query_string = '(ITS1 OR ITS2 OR (internal transcribed spacer) OR 18S OR SSU OR LSU OR 28S OR 5.8S) AND Fungi[Organism]';
my $gb = Bio::DB::GenBank->new();
mkdir($basedir) unless -d $basedir;

my $query;
for my $year ( $mindate[0]..$maxdate[0]) {
   my $ofile = File::Spec->catfile($basedir,"$year\_ITS.gbk");
   next if( -f $ofile && ! -z $ofile);
   my $query = Bio::DB::Query::GenBank->new(-db=>'nucleotide',
                                         -query =>$query_string,
					 -mindate=>join("/",$year,1,1),
					 -maxdate=> join("/",$year,12,31),
					-verbose => $DEBUG,
                                         );
    print "Total records for $year query is ", $query->count,"\n";
    my $out = Bio::SeqIO->new(-file => ">$ofile",-format=>'genbank');
    my $stream = $gb->get_Stream_by_query($query);
    while (my $seq = $stream->next_seq) {
     $out->write_seq($seq);
   }
}
