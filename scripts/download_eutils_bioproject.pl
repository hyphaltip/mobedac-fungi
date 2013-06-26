#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use LWP::Simple;
use Encode;
use Cache::File;
use XML::Simple;
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Env qw(USER);

my $SLEEP_TIME = 2;
my $cache_dir = "/tmp/eutils_".$ENV{USER};
my $cache_filehandle;
my $basedir = 'genomes_download';

my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
my $query='txid4751[Organism:exp]';

my $force = 0;
my $debug = 0;
my $retmax = 1000;

my $use_cache = 1;
GetOptions(
    'q|query:s' => \$query,
    'b|basedir:s' => \$basedir,
    'debug|v!' => \$debug,
    'retmax:i'  => \$retmax,
    'f|force!'  => \$force, # force downloads even if file exists
    'cache!'    => \$use_cache,
    );

$SLEEP_TIME = 0 if $debug; # let's not wait when we are debugging

if( $use_cache ) {
    &init_cache();
}

my $db = 'genome';
my $url = sprintf("esearch.fcgi?db=%s&term=%s&rettype=acc&retmax=%d&usehistory=y",
		  $db,$query,$retmax);

mkdir($basedir) unless -d $basedir;

warn "url is $url\n" if $debug;
my $output = get_web_cached($base,$url);
my ($web,$key,$count);

if( $output =~ /<WebEnv>(\S+)<\/WebEnv>/) {
    $web = $1;
}
if ($output =~ /<QueryKey>(\d+)<\/QueryKey>/) {
    $key = $1;
}
if( $output =~ /<Count>(\d+)<\/Count>/) {
    $count = $1;
}

my @ids;
while ($output =~ /<Id>(\d+?)<\/Id>/sg) {
   push(@ids, $1);
}

warn "processing ", scalar @ids," project IDs\n";

for my $id ( @ids ) {    
    warn( "id is $id\n") if $debug;
    my $xs = XML::Simple->new;

    $url = sprintf('esummary.fcgi?db=genome&id=%d',$id);    
    # post the elink URL
    $output = get_web_cached($base,$url);
    my $simplesum = $xs->XMLin($output);
    my $doc = $simplesum->{DocSum};
    my $species;
    my @keyorder;

    for my $item ( @{$doc->{Item} } ) {
	next unless defined $item->{content} ;
	push @keyorder,[ $item->{Name}, $item->{content} ];
	if( $item->{Name} eq 'Organism_Name' ) {
	    $species = $item->{content};
	    $species =~ s/\s+/_/g;
	}
    }
    if( ! defined $species ) {
	warn("Cannot parse species from $id\n");
	$species = 'Unknown';
    }
    my $species_dir = File::Spec->catdir($basedir,$species);
    mkdir($species_dir);
    print "working species $species\n";
    open(my $summary_rpt => ">$species_dir/SPECIES_REPORT.txt") || die "Cannot open $species_dir/SPECIES_REPORT.txt";
    for my $key ( @keyorder ) {
	print $summary_rpt join("\t", @$key), "\n";
    }
    close($summary_rpt);

    $url = sprintf('elink.fcgi?dbfrom=genome&db=bioproject&id=%d&retmax=1000',$id);
    warn( "id is $id\n") if $debug;
    # post the elink URL
    $output = get_web_cached($base,$url);
    my %seen_projids;
    if( -f "$species_dir/bioproject_ids.dat" ) {
	open(my $fh => "$species_dir/bioproject_ids.dat");
	while(<$fh>) {
	    next if /^BIOPROJECT/;
	    my ($pid,$genomic) = split;
	    $seen_projids{$pid} = $genomic;
	}
    } 
    my @bioproj_ids;
    if( $output =~ /<LinkSetDb>(.+)<\/LinkSetDb>/sio ) {
	my $link = $1;
	while ($link =~ /<Id>(\d+?)<\/Id>/sg) {
	    if( ! exists $seen_projids{$1} ) {
		$seen_projids{$1} = -1;
	    }
	}
	my %info;
	for my $bioproj_id ( sort { $a <=> $b } keys %seen_projids ) {	    
	    warn("projid is $bioproj_id and seen status is $seen_projids{$bioproj_id}\n") if $debug;
	    
	    next unless ( $force || $seen_projids{$bioproj_id} < 0 );
	    $url = sprintf('esummary.fcgi?db=bioproject&id=%d',$bioproj_id);
	    $output = get_web_cached($base,$url);
	    warn("obtained from $url\n") if $debug;
	    warn("$output") if $debug;

	    my $simple = $xs->XMLin($output,ForceArray => ['Project_Objectives_List']);
	    my ($projname, $projtitle,$data_type,$objective);
	    if( ! ref($simple) ) { 
		warn("could not parse summary doc\n");
		next;
	    }
	    # assume a single document summary will be present since we are
	    # requesting this fetch for a single bioproject
	    my $doc = $simple->{DocumentSummarySet}->{DocumentSummary};
	    
	    unless( $data_type = $doc->{Project_Data_Type} ) {
		die("cannot data type description out of the XML\n");
	    }

	    warn("data type is $data_type\n") if $debug;

	    OBJ: for my $obj ( @{$doc->{Project_Objectives_List}} ) {
		my $struct = $obj->{Project_Objectives_Struct};
		if( ref($struct) =~ /ARRAY/ ) {
		    for my $n ( @{$struct} ) {
			my $r = $n->{Project_ObjectivesType};
			$objective = $r;
			last OBJ if( $r =~ /Assembly|Annotation|RefSeq/ );
		    }
		} elsif( ref($struct) =~ /HASH/ ) {
		    my $r = $struct->{Project_ObjectivesType};
		    $objective = $r;
		    last OBJ if( $r =~ /Assembly|Annotation|RefSeq/ );
		}
	    }
	    unless( $projtitle = $doc->{Project_Title} ) {
		warn("cannot parse description out of the XML\n");
	    } 
	    unless( $projname = $doc->{Organism_Name} ) {
		# OrganismName
		# Strain
		warn("cannot parse name out of the XML\n");
	    }
	    warn("objective is $objective\n") if $debug;
	    if( $objective =~ /Assembly|Annotation|RefSeq/ ||
		$data_type =~ /RefSeq/) {
		$seen_projids{$bioproj_id} = 1;
	    } else {
		$seen_projids{$bioproj_id} = 0;
		next;
	    }
	    
	    $projname =~ s/[\s+\/\\'"]/_/g;
	    push @{$info{$projname}},    {id        => $bioproj_id,
					  objective => $objective,
					  name      => $projname,
					  type      => $data_type,
					  title     => $projtitle};	    
	}
	
	for my $name ( keys %info ) {
	    my @process;
	    my $saw_ref = 0;
	    for my $inf ( @{$info{$name}} ) {
		if( $inf->{type} =~ /mitochondria/i ) {
		    warn("mito is $name\n");
		    push @process, $inf;
		} elsif( $inf->{type} =~ /RefSeq|Reference/i ) {
		    $saw_ref = 1;
		    push @process, $inf;			
		}
	    }
	    if( ! $saw_ref ) {
		# take it all - the only issue is if there are multiple assemblies for same strain, 
		# but this will require human intervention
		@process = @{$info{$name}};	    
	    }	

	    for my $inf ( @process ) {
		my $is_ref = ( $inf->{type} =~ /RefSeq|Reference/i );
		$url = sprintf('elink.fcgi?dbfrom=bioproject&db=nuccore&id=%d&retmax=1000',
			       $inf->{id});
		warn("Getting ",$inf->{id}," $species $url\n") if $debug;
		
		my $output2 = get_web_cached($base,$url);
		warn("obtained from $url\n$output2\n") if $debug;
		# parse WebEnv and QueryKey
		my @nucl_ids;
		my $links;
		if( $output2 =~ /<LinkSetDb>(.+)<\/LinkSetDb>/si) {
		    $links = $1;	
		    while($links =~ /<Id>(\d+?)<\/Id>/sg ) {
			push @nucl_ids, $1;
		    }
		    $seen_projids{$inf->{id}} = scalar @nucl_ids;
		    warn("links were $links (@nucl_ids\n") if $debug;
		    next unless( @nucl_ids );
		    
		    my $species_strain_dir = File::Spec->catdir($species_dir,$inf->{name});
		    # convert this into creating XML for EUPATHDB?
		    if( -d $species_strain_dir && ! $force ) {
			warn("skipping $species_strain_dir, already processed\n");
			next;
		    }
		    mkdir($species_strain_dir);

		    my $rptfile = File::Spec->catfile($species_strain_dir,$inf->{id}.".txt");
		    open(my $rptfh => ">$rptfile") || die "$rptfile: $!";
		    print $rptfh join("\n",$inf->{title}, $inf->{name}, $inf->{objective},
				      sprintf("Is%sreference",$is_ref ? ' ' : ' not ')),"\n";
		    
		    warn Dumper($doc->{Lineage}) if $debug;
		    my (@l,@ln);
		    for my $item ( @{$doc->{Lineage}->{Item} || []} ) {
			push @l, $item->{taxid};
			push @ln, $item->{content};
#		print $item->{taxid}, " ", $item->{content}, "\n";
		    }
		    print $rptfh join(";", @l), "\n";
		    print $rptfh join(";", @ln), "\n";
		    
		    print $rptfh join("\t", qw(GENOMEGROUP BIOPROJECT IS_REF NUCL_GI_IDS)),"\n";		
		    print $rptfh join("\t", $id, $inf->{id}, join(",", @nucl_ids)),"\n";
		}
	    }
	}
	open(my $projidsfh => ">$species_dir/bioproject_ids.dat") || die "cannot open proj_id file: $!";
	print $projidsfh join("\t", qw(BIOPROJECT_ID IS_GENOMIC)),"\n";
	if( keys %seen_projids ) {
	    print $projidsfh join("\n", map { join("\t", $_, $seen_projids{$_}) } 
				  sort { $a <=> $b }  
				  keys %seen_projids), "\n";
	    
	}
	close($projidsfh);
	last if $debug;
    }
    last if $debug;
}

sub init_cache {
    if( ! $cache_filehandle ) {
	mkdir($cache_dir) unless -d $cache_dir;	
	$cache_filehandle = Cache::File->new( cache_root => $cache_dir);
    }
}

sub get_web_cached {
    my ($base,$url) = @_;
    if( ! defined $base || ! defined $url ) {
	die("need both the URL base and the URL stem to proceed\n");
    }
    unless( $use_cache ) {
	sleep $SLEEP_TIME;
	return get($base.$url);
    }
    my $val = $cache_filehandle->get($url);
    unless( $val ) {
	$val = encode("utf8",get($base.$url));
	sleep $SLEEP_TIME;
	$cache_filehandle->set($url,$val,'1 day');
    } 
    return decode("utf8",$val);    
}
