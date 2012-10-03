#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use Bio::Tools::EUtilities::History;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;
use File::Temp qw(tempfile);

# THIS REQUIRES YOU SHOULD HAVE RUN download_eutils_bioproject.pl already
my $MAX_TO_QUERY = 500;
my $basedir = 'genomes_download';

my $force = 0;
my $debug = 0;
my $retmax = 1000;
my $skip_mRNA = 1;
my $tmpdir = '/dev/shm';
GetOptions(
    'b|basedir:s' => \$basedir,
    'debug|v!' => \$debug,
    'retmax:i'  => \$retmax,
    's|skip|skipmRNA!' => \$skip_mRNA,
    'f|force!'  => \$force, # force downloads even if file exists	
    't|tmp:s'  => \$tmpdir,
    );

my ($tempfh,$tempfile) = tempfile(DIR => $tmpdir);
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

	    warn("downloading @seqids\n") if $debug;		
	    my $first = 1;
	    my $skip_project = 0;
	    if( @seqids ) {
		warn("processing $dir/$projectid\n");
	    }
            my @needed_ids;
	    for my $id ( @seqids ) {
		if( $skip_project ) {
		    warn("Skipping out on project $projectid, it is mRNA\n");
		    last;
		}
		my $outfile = "$basedir/$dir/$projectid/$id.gbk";
		if( $force || ! -f $outfile ) {
		    push @needed_ids, $id;
	 	}
		if( -f $outfile ) {
		    if( $first && $skip_mRNA ) {
			open(my $fh => $outfile)|| die $!;
			while(<$fh>) {
			    if( /^LOCUS/ && /mRNA/ ) {
				$skip_project = 1;
				last;
			    }
			}
			close($fh);
		    }
		    $first = 0;
		}
	    }
	    next if( $skip_project );
	    while( @needed_ids ) {
		my @chunk = splice(@needed_ids,0,$MAX_TO_QUERY);
		warn("chunk is @chunk\n") if $debug;
		my %chunkids = map { $_ => 0 } @chunk;
		my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
						       -db      => 'nucleotide',
						       -rettype => 'gbwithparts',
						       -email   => 'jason@bioperl.org',
						       -id      => \@chunk);		    
#		$factory->set_parameters(-rettype => 'gbwithparts');
		eval {
		    $factory->get_Response(-file => $tempfile);
		    
		    my $in = Bio::SeqIO->new(-format => 'genbank', 
					     -file => $tempfile);
		    $first = undef;
		    while( my $seq = $in->next_seq ) {
			my $gi = $seq->primary_id;
			$chunkids{$gi}++;
			my $out = Bio::SeqIO->new
			    (-format => 'genbank', 
			     -file =>">$basedir/$dir/$projectid/$gi.gbk");
			$first = "$basedir/$dir/$projectid/$gi.gbk" unless defined $first;
			$out->write_seq($seq);
		    }
		    $in->close();
		    open(my $fh => $first)|| die $!;
		    while(<$fh>) {
			if( /^LOCUS/ && /mRNA/ ) {
			    $skip_project = 1;
			    last;
			}
		    }
		    for my $seqid ( keys %chunkids ) {
			if( $chunkids{$seqid} == 0 ) {
			    warn("Download proble with $dir/$projectid/$seqid\n");
			}
		    }
		};
		if( $@ ) {
		    warn("Download problem with $dir/$projectid from: $@\n");
		}
		if( $skip_project ) {
		    warn("Skipping out on project $projectid, it is mRNA\n");
		    last;
		}
		sleep 3;
	    }
	}
    }
    #last if $debug;
}
