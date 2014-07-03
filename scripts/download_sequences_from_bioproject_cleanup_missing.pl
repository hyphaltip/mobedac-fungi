#!/usr/bin/perl -w
use strict;
use Bio::DB::EUtilities;
use Bio::SeqIO;
use Getopt::Long;
use File::Copy qw(copy);
use File::Temp qw(tempfile);
# let's find the 
my $MAX_TO_QUERY = 100;
my $basedir = 'fungal_genomes_download';

my $force = 0;
my $debug = 0;
my $runone = undef;
my $retmax = 1000;
my $skip_mRNA = 1;
my $tmpdir = '/dev/shm';
my $match;
GetOptions(
    'b|basedir:s' => \$basedir,
    'debug|v!' => \$debug,
    'runone:s'  => \$runone,
    'retmax:i'  => \$retmax,
    's|skip|skipmRNA!' => \$skip_mRNA,
    'f|force!'  => \$force, # force downloads even if file exists	
    't|tmp:s'  => \$tmpdir,
    'max:i'    => \$MAX_TO_QUERY,
    'match:s'  => \$match,
    );

my ($tempfh,$tempfile) = tempfile(DIR => $tmpdir,CLEANUP=>1);
opendir(BASE,$basedir) || die "cannot open $basedir directory: $!";
for my $dir ( readdir(BASE) ) {
    next if $dir =~ /^\./;
    next unless ( -d "$basedir/$dir");
    if( $match && $dir !~ /$match/ ) {
	next;
    }
    if( $runone ) {
	next unless $dir =~ /$runone/;
    }
    # each dir is a species name
    opendir(DIR,"$basedir/$dir") || die "cannot open $basedir/$dir: $!";
    for my $projectid ( readdir(DIR) ) {
	next if $projectid =~ /^\./;
	next unless ( -d "$basedir/$dir/$projectid");	
	opendir(PROJDIR,"$basedir/$dir/$projectid") || die "cannot open $basedir/$dir/$projectid: $!";
	
	my %filenames;
	my $species_name_base;
	my $skip_project = 0;
	for my $pfile ( readdir(PROJDIR)  ) {
	    if(  $pfile =~ /^(\S+)\.txt/ ) {		
		my $nm = $1;
		open(my $fh => "$basedir/$dir/$projectid/$pfile") || die $!;
		my $projtitle = <$fh>;
		my $project_id_str = <$fh>;
		my $objective = <$fh>;
		my $is_reference = <$fh>;
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
		close($fh);
		$species_name_base = $nm;
	    } else {
		if( ! -z "$basedir/$dir/$projectid/$pfile" ) {
		    $filenames{$pfile} = "$basedir/$dir/$projectid/$pfile";
		}
	    }
	}
	unless( $species_name_base ) {
	    warn("no project file with species name ready to process in $dir/$projectid\n");
	    next;
	}
	my $master_file = "$basedir/$dir/$projectid/$species_name_base" . '.master.genbank';
	# master exists, skip
	next unless ( $force || ! $filenames{$master_file} );

	my (@acc, @acc_scaff, %skip_toplevel);
	my $molecule = 'nuclear';
	# %skip_toplevel - is just to remember to skip this file in the future as it is a toplevel 
	# organizing sequence file
	my %mt;
	for my $gbk ( map { $_->[1] } 
		      sort { $a->[0] <=> $b->[0] }
		      map { [ /(\d+)\.gbk/, $_ ] } 
		      grep { /\d+\.gbk$/ } keys %filenames ) {
	    open(my $infh => $filenames{$gbk} ) || die $!;
	    while(<$infh>) {
		if( $skip_mRNA &&
		    /^LOCUS/ && /mRNA/ ) {
		    $skip_project = 1;
		    last;
		}
		if(/^WGS\s+(\S+)/) {
		    my $accnums = $1;
		    push @acc,split(/,/,$accnums);
		    $skip_toplevel{$gbk}++;
		} elsif( /^WGS_SCAFLD\s+(\S+)/ ) {
		    my $accnums = $1;
		    push @acc_scaff,split(/,/,$accnums);
		    $skip_toplevel{$gbk}++;
		} elsif( /organelle="mitochondrion"/ ) {
		    $mt{$gbk} = $filenames{$gbk};		    
		} elsif( /ORIGIN/ ) {
		    last;
		}
	    }
	    last if $skip_project;
	}
	if( $skip_project ) {
	    warn("skipping mRNA project $dir -> $projectid\n");
	    next;
	}
	
	my @set = scalar @acc_scaff > 0 ? @acc_scaff : @acc;
	
	if( @set ) {
	    warn( "$dir $projectid ACC are ", join("\n", @set),"\n");	
	    my %full_acc_list;
	    for my $acc ( @set ) {
		if( $acc =~ /(\S+)\-(\S+)/ ) {
		    warn("$acc\n") if $debug;
		    my ($start,$end) = ($1,$2);
		    my $break_point = 0;
		    for( my $i =0; $i < length($start); $i++) {
			if( substr($end,$i,1) ne substr($start,$i,1) ) {
			    $break_point = $i;
			    last;
			}
		    }
		    my $start_num = substr($start,$break_point-1);
		    my $end_num   = substr($end,$break_point-1);
		    my $l = length $start_num;
		    my @acc_full;
		    warn("$start $start_num -> $end_num\n") if $debug;
		    for(my $j = $start_num; $j <= $end_num; $j++ ) {			
			$full_acc_list{sprintf "%s%0$l"."d",
				       substr($start,0,$break_point-1),
				       $j}++;
		    }
		} else {
		    $full_acc_list{$acc}++;
		}
	    }
	    my @list_ids =sort keys %full_acc_list;
	    while( @list_ids ) {
		my @chunk = splice(@list_ids, 0, $MAX_TO_QUERY);
		warn("downloading @chunk\n") if $debug;
		my $factory = Bio::DB::EUtilities->new(-eutil   => 'efetch',
						       -db      => 'nucleotide',
						       -id      => \@chunk,
						       -email   => 'jason@bioperl.org',
						       -rettype => 'gbwithparts');
		eval { 
		    $factory->get_Response(-file => ">$tempfile");		    
		    local($/) = "\/\/\n";
		    #my $in = Bio::SeqIO->new(-format => 'genbank', -file => $tempfile);
		    open(my $infh => $tempfile) || die $!;
		    while(<$infh>) {
			next unless length ($_) > 0;
			next if /^\s+$/;
			my $accn = '?';		
			if( /^ACCESSION\s+(\S+)/ || /LOCUS\s+(\S+)/) {
			    $accn = $1;
			}
			if( /GI:(\d+)/ ) {
			    my $ginum = $1;	
			    if( ! -f "$basedir/$dir/$projectid/$ginum.gbk" ) {
				open(my $tfh => ">$basedir/$dir/$projectid/$ginum.gbk") || die $!;
				print $tfh $_;
				close($tfh);
				$filenames{"$ginum.bak"} = "$basedir/$dir/$projectid/$ginum.gbk";
			    } else {
				warn("not writing $accn -> $ginum, already exists\n") if $debug;
			    }
			} else { 
			    warn("no GI number parseable ($_)\n");
			}
		    }
#		while(my $seq = $in->next_seq ) {
#		my $ginum = $seq->primary_id;
#		if( ! -f "$basedir/$dir/$projectid/$ginum.gbk" ) {
#		    Bio::SeqIO->new(-format => 'genbank',
#				    -file   => ">$basedir/$dir/$projectid/$ginum.gbk")->write_seq($seq);
#		    $filenames{"$ginum.bak"} = "$basedir/$dir/$projectid/$ginum.gbk";
#		} else {
#		    warn("not writing ", $seq->accession_number, " -> $ginum, already exists\n");
#		}
		};			    
		if( $@ ) {
		    warn("@chunk, $@\n");
		}		
	    }
	}
	open(my $master_out => ">$master_file") || die $!;
	for my $gbk ( map { $_->[1] } 
		      sort { $a->[0] <=> $b->[0] }
		      map { [ /(\d+)\.gbk/, $_ ] } 
		      grep { /\d+\.gbk$/ } keys %filenames ) {
	    next if $skip_toplevel{$gbk} || $mt{$gbk};
	    open(my $infh => $filenames{$gbk} ) || die $!;
	    while(<$infh>) {
		print $master_out $_;
	    }
	    close($infh);
	}
	open(my $master_mt => ">$basedir/$dir/$projectid/$species_name_base" . '.MT.genbank') || die $!;
	for my $mt ( keys %mt ) {
	    open(my $infh => $mt{$mt}) ||die $!;
	    while(<$infh>) {
		print $master_mt $_;
	    }
	    close($infh);
	}
    }
    last if $runone;
}
