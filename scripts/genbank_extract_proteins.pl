#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use Bio::SeqIO;

my $debug;
GetOptions(
    'v|verbose|debug!' => \$debug,
    );

my $dir = shift || die $!;

opendir(DIR,$dir) || die $!;

for my $file ( readdir(DIR) ) {
    next unless $file =~ /(\S+)\.genbank\.gz$/;
    my $out = Bio::SeqIO->new(-format => 'fasta',
			      -file   => ">$dir/$1.peps.fa");
    open(my $fh => "zcat $dir/$file |") || die $!;
    my $in = Bio::SeqIO->new(-format => 'genbank',
			     -fh   => $fh) || die $!;
    while(my $seq = $in->next_seq ) {
	for my $cds ( grep { $_->primary_tag eq 'CDS' } 
		      $seq->get_SeqFeatures ) 
	{	    
	    my ($name,$gi,$acc,$product);
	   
	    if( $cds->has_tag('protein_id') ) {
		($acc) = $cds->get_tag_values('protein_id');
	    } else { warn ("no protein_id $name ", $seq->display_id, " ", $file,"\n"); 
		     $acc = '';
	    }
	    if( $cds->has_tag('db_xref') ) {
		($gi) = $cds->get_tag_values('db_xref');
	    } else { warn ("no db_xref $name $acc", $seq->display_id, " ", $file,"\n"); 
		     $gi = '';
	    }	    
	    if( $cds->has_tag('locus_tag') ) {
		($name) = $cds->get_tag_values('locus_tag');
	    } else { 
		if( $cds->has_tag('gene'))  {
		    ($name) = $cds->get_tag_values('gene');
		} else {
		    warn("no locus_tag or gene ", $seq->display_id, " ", $file,"\n");
		    $name = '';
		}
	    }
	    $name ||= $acc || $gi;
	    my $desc = sprintf("gi=%s acc=%s loc=%s:%s",
			       $gi,$acc,$seq->display_id,
			       $cds->location->to_FTstring);
	    if( $cds->has_tag('product') ) {
		$desc .= sprintf(" product=\"%s\"",
				 join("; ",$cds->get_tag_values('product')));
	    }
	    if( $cds->has_tag('note') ){
		$desc .= sprintf(" note=\"%s\"",
				 join("; ", $cds->get_tag_values('note')));
	    }
	    if( $cds->has_tag('translation') ) {
		my ($pepseq) = $cds->get_tag_values('translation');
		my $pep = Bio::Seq->new(-seq => $pepseq,
					-id => $name,
					-desc =>  $desc);
		$out->write_seq($pep);
	    } else {
		my ($codon_start) = 1;
		if( $cds->has_tag('codon_start') ) {
		    ($codon_start) = $cds->get_tag_values('codon_start');
		}
		my $transl_table = 1;
		if( $cds->has_tag('transl_table') ) {
		    ($transl_table) = $cds->get_tag_values('transl_table');
		}
		my $pep = $cds->spliced_seq->translate
		    (-frame => $codon_start -1,
		     -codontable_id => $transl_table);
		$pep->id($name);
		
		$pep->desc($desc);
		$out->write_seq($pep);
	    }
			
	    last if $debug;
	}
	last if $debug;
    }
#    last if $debug;
}
