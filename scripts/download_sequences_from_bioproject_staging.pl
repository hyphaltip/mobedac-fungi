#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use File::Spec;
use Getopt::Long;
use File::Temp qw(tempfile);

# THIS REQUIRES YOU SHOULD HAVE RUN download_eutils_bioproject.pl already
my $MAX_TO_QUERY = 500;
my $basedir = 'fungal_genomes_download';

my $force = 0;
my $debug = 0;
my $retmax = 1000;
my $runone = 0;
my $skip_mRNA = 1;
my $tmpdir = '/dev/shm';
GetOptions(
    'b|basedir:s' => \$basedir,
    'debug|v!' => \$debug,
    'retmax:i'  => \$retmax,
    's|skip|skipmRNA!' => \$skip_mRNA,
    'f|force!'  => \$force, # force downloads even if file exists	
    't|tmp:s'  => \$tmpdir,
    'runone!'  => \$runone,
    );

my ($tempfh,$tempfile) = tempfile(DIR => $tmpdir);
opendir(BASE,$basedir) || die "cannot open $basedir directory: $!";
for my $species_dir ( readdir(BASE) ) {
    next if $species_dir =~ /^\./;
    next unless ( -d "$basedir/$species_dir");

    # each dir is a species name
    opendir(DIR,"$basedir/$species_dir") || die "cannot open $basedir/$species_dir: $!";
    for my $strain_dir ( readdir(DIR) ) {
	next if $strain_dir =~ /^\./;
	next unless ( -d "$basedir/$species_dir/$strain_dir");
	opendir(PROJDIR,"$basedir/$species_dir/$strain_dir") || die "cannot open $basedir/$species_dir/$strain_dir: $!";
	for my $projfile ( readdir(PROJDIR)  ) {
	    next unless $projfile =~ /^(\S+)\.txt/;
	    my $nm = $1;
	    open(my $fh => "$basedir/$species_dir/$strain_dir/$projfile") || die $!;
	    my $projtitle = <$fh>;
	    my $project_id_str = <$fh>;	    
	    my $objective = <$fh>;
	    my $is_reference = <$fh>;
	    my $header;
	    while(<$fh>) {
		$header = $_;
		last if $header =~ /^GENOMEGROUP/;
	    }
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
		warn("processing $species_dir/$strain_dir $projfile\n");
	    }
            my @needed_ids;
	    for my $id ( @seqids ) {
		if( $skip_project ) {
		    warn("Skipping out on project $species_dir/strain_dir, it is mRNA\n");
		    last;
		}
		my $outfile = "$basedir/$species_dir/$strain_dir/$id.gbk";
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
			     -file =>">$basedir/$species_dir/$strain_dir/$gi.gbk");
			$first = "$basedir/$species_dir/$strain_dir/$gi.gbk" unless defined $first;
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
			    warn("Download proble with $species_dir/$strain_dir/$seqid\n");
			}
		    }
		};
		if( $@ ) {
		    warn("Download problem with $species_dir/$strain_dir from: $@\n");
		}
		if( $skip_project ) {
		    warn("Skipping out on project $strain_dir, it is mRNA\n");
		    last;
		}
		sleep 3;
	    }
	}
    }
    last if $runone;
}
