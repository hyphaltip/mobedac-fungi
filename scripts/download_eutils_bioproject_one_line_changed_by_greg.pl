#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use Bio::Tools::EUtilities::History;
use LWP::Simple;
use XML::Simple;
use File::Spec;
use Getopt::Long;
use Data::Dumper;
use Bio::DB::Taxonomy;
use Bio::DB::GenBank;

my $basedir = 'genomes_download';

my $taxdb = Bio::DB::Taxonomy->new(-source => 'entrez');
#my $xs = XML::Simple->new;
my $gb = Bio::DB::GenBank->new;

my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils';
my $query='txid4751[Organism:exp]';

my $force = 0;
my $debug = 0;
my $retmax = 1000;

GetOptions(
    'q|query:s' => \$query,
    'b|basedir:s' => \$basedir,
    'debug|v!' => \$debug,
    'retmax:i'  => \$retmax,
    'f|force!'  => \$force, # force downloads even if file exists
    );


my $db = 'genome';
my $url = sprintf("%s/esearch.fcgi?db=%s&term=%s&rettype=acc&retmax=%d&usehistory=y",
		  $base, $db,$query,$retmax);

mkdir($basedir) unless -d $basedir;

warn "url is $url\n" if $debug;
my $output = get($url);
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
#print "$#ids\n";

#warn "processing ", scalar @ids," project IDs\n";
my $size=@ids;
warn "processing ", $size," project IDs\n";
#print "$output";
#exit;

for my $id ( @ids ) {    
#next if $id ne "9518";
    warn( "id is $id\n") if $debug;

    my $url = sprintf('%s/esummary.fcgi?db=genome&id=%d',$base,$id);    
    # post the elink URL
    $output = get($url);
    my $xs = XML::Simple->new;
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

    $url = sprintf('%s/elink.fcgi?dbfrom=genome&db=bioproject&id=%d&retmax=1000',$base,$id);
    warn( "id is $id\n") if $debug;
    # post the elink URL
    $output = get($url);
    my %seen_projids;
    if( -f "$species_dir/bioproject_ids.dat" ) {
     open(my $fh => "$species_dir/bioproject_ids.dat");
     while(<$fh>) {
	my ($pid) = split;
	$seen_projids{$pid}++;
    }
   } 
    my @bioproj_ids;
    if( $output =~ /<LinkSetDb>(.+)<\/LinkSetDb>/sio ) {
	my $link = $1;
	while ($link =~ /<Id>(\d+?)<\/Id>/sg) {
	    push(@bioproj_ids, $1);
	}
        open(my $projidsfh => ">$species_dir/bioproject_ids.dat") || die $!;
	print $projidsfh join("\n", sort { $a <=> $b } @bioproj_ids), "\n";
	close($projidsfh);
	for my $bioproj_id ( @bioproj_ids ) {
	    next unless ( $force || ! exists $seen_projids{$bioproj_id} || ! $seen_projids{$bioproj_id});
	    $url = sprintf('%s/efetch.fcgi?db=bioproject&id=%d',$base,$bioproj_id);
	    $output = get($url);
	    $url = sprintf('%s/elink.fcgi?dbfrom=bioproject&db=nuccore&id=%d&retmax=1000',
			   $base,$bioproj_id);
	    warn("Getting $bioproj_id $species $url\n") if $debug;
	    my $output2 = get($url);
	    warn("obtained from $url\n") if $debug;
	    #print $output2;
	    #next;
	    #warn("$output");
	    # parse WebEnv and QueryKey
	    my @nucl_ids;
	    my $links;
	    if( $output2 =~ /<LinkSetDb>(.+)<\/LinkSetDb>/si) {
		$links = $1;	
		while($links =~ /<Id>(\d+?)<\/Id>/sg ) {
		    push @nucl_ids, $1;
		}
		next unless @nucl_ids;
		
		my $simple = $xs->XMLin($output);
		my ($projname, $projtitle);
		if( ! ref($simple) ) { 
		   next;
		}
		my $doc = $simple->{DocumentSummary};
		
		unless( $projtitle = $doc->{Project}->{ProjectDescr}->{Title} ) {
		    warn("cannot parse name out of the XML\n");
		} 
		unless( $projname = $doc->{Project}->{ProjectDescr}->{Name} ) 
		{	# OrganismName
		    # Strain
		    warn("cannot parse name out of the XML\n");
		}
		my $species_strain_dir = File::Spec->catdir($species_dir,$bioproj_id);
		# convert this into creating XML for EUPATHDB?
		if( -d $species_strain_dir && ! $force ) {
			warn("skipping $species_strain_dir, already processed\n");
			next;
		}
		mkdir($species_strain_dir);
		#$projname =~ s/\s+/_/g;
		$projname =~ s/(\s+|\/)/_/g; #for some case the name of the project contain /
		open(my $rptfh => ">$species_strain_dir/$projname.txt") || die "$species_strain_dir/$projname.txt: $!";
		print $rptfh $projtitle, "\n", $projname, "\n";	    
		
		warn Dumper($doc->{Lineage}) if $debug;
		my (@l,@ln);
		for my $item ( @{$doc->{Lineage}->{Item} || []} ) {
		    push @l, $item->{taxid};
		    push @ln, $item->{content};
#		print $item->{taxid}, " ", $item->{content}, "\n";
		}
		print $rptfh join(";", @l), "\n";
		print $rptfh join(";", @ln), "\n";
		
		print $rptfh join("\t", qw(GENOMEGROUP BIOPROJECT NUCL_GI_IDS)),"\n";		
		print $rptfh join("\t", $id, $bioproj_id, join(",", @nucl_ids)),"\n";
		if( 0 ) {
		    if( $force || 
			! -f "$species_strain_dir/sequences.gbk" ||
			-z "$species_strain_dir/sequences.gbk" ) {
			eval { 
			    my $out = Bio::SeqIO->new(-format => 'genbank',
						      -file   => ">$species_strain_dir/sequences.gbk");
			    for my $nuclid ( @nucl_ids ) {
				if( my $seq = $gb->get_Seq_by_gi($nuclid) ) {
				    $out->write_seq($seq);
			    } else {
				print $rptfh "Seq $nuclid failed\n";
			    }
			    }
			};
			if( $@ ) {
			warn("error $@");
			}
		    } else {
			warn("skipping $species_strain_dir sequences, already downloaded\n");
		    }
		}
#		sleep 5;
	    }
	} 
#	sleep 5;
    }
#    sleep 10;
}
