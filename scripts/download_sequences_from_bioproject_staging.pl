#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use Bio::Tools::EUtilities::History;
use Bio::DB::GenBank;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Bio::DB::Query::GenBank;

# THIS REQUIRES YOU SHOULD HAVE RUN download_eutils_bioproject.pl already

my $basedir = 'genomes_download';

my $gb = Bio::DB::GenBank->new();

my $force = 0;
my $debug = 0;
my $retmax = 1000;

GetOptions(
    'b|basedir:s' => \$basedir,
    'debug|v!' => \$debug,
    'retmax:i'  => \$retmax,
    'f|force!'  => \$force, # force downloads even if file exists
    );



opendir(BASE,$basedir) || die "cannot open $basedir directory: $!";
for my $dir ( readdir(BASE) ) {
    next if $dir =~ /^\./;
    next unless ( -d "$basedir/$dir");
    
    # each dir is a species name
    opendir(DIR,"$basedir/$dir") || die "cannot open $basedir/$dir: $!";
    for my $projectid ( readdir(DIR) ) {
	next if $projectid =~ /^\./;
	next unless ( -d "$basedir/$dir/$projectid");
	opendir(PROJFILE,"$basedir/$dir/$projectid") || die "cannot open $basedir/$dir/$projectid: $!";
	for my $projfile ( readdir(PROJFILE)  ) {
	    next unless $projfile =~ /^(\S+)\.txt/;
	    my $nm = $1;
	    open(my $fh => "$basedir/$dir/$projectid/$projfile") || die $!;
	    my $projtitle = <$fh>;
	    my $project_id_str = <$fh>;
	    my $taxonomyid_str = <$fh>;
	    my $taxonomy_str = <$fh>;
	    my $header = <$fh>;
	    chomp($projectid);
	    my @header = split(/\s+/,$header);
	    my @seqids;
	    while(<$fh>) {
		my ($genomegroup,$bioproject,$nuclids) = split;
		push @seqids, split(/,/,$nuclids);
	
	    }

	    warn("downloading @seqids\n");		
	    for my $id ( @seqids ) {
		my $outfile = "$basedir/$dir/$projectid/$id.gbk";
		if( $force || ! -f $outfile ) {		
		    my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
							   -db      => 'nucleotide',
							   -rettype => 'gb',
							   -email   => 'jason.stajich@gmail.com',
							   -id      => $id);		    
		    $factory->set_parameters(-rettype => 'gbwithparts');
		    eval {
			$factory->get_Response(-file => $outfile);
		    };

		    if( $@ ) {
			warn("Download problem with $id from: $@\n");
		    }
		}	
	    } else {
		warn("already see $basedir/$dir/$projectid/$nm.gbk present so skipping\n");
	    } 
	}
    }
    last if $debug;
}
