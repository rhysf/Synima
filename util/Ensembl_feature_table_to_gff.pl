#!/usr/bin/perl -w
use strict;
use warnings;

# Usage
my $file = shift or die "Usage: $0 <file> > genome.annotation.gff3\n";
open (FILE, "<$file") or die "Failed to open file\n";

while (<FILE>) {
    chomp;
    unless ( m/^\#/) { 
	my ($feature,
	    $class,
	    $assembly,
	    $assembly_unit,
	    $seq_type,
	    $chromosome,
	    $genomic_accession,
	    $start,
	    $end,
	    $strand,
	    $product_accession,
	    $non_redundant_refseq,
	    $related_accession,
	    $name,
	    $symbol,
	    $GeneID,
	    $locus_tag,
	    $feature_interval_length,
	    $product_length,
	    $attributes
	    )= split /\t/;


	my $score='.';
	my $phase = '.';
	my $source = '.';

	if ($feature eq 'CDS' and length $product_accession ) {
	    
	    print "$genomic_accession";
	    print "\t";
	    print "$source";
	    print "\t";
	    print "mRNA";
	    print "\t";
	    print "$start";
	    print "\t";
	    print "$end";
	    print "\t";
	    print "$score";
	    print "\t";
	    print "$strand";
	    print "\t";
	    print "$phase";
	    print "\t";
	    print "ID=$product_accession";
	    print "\n";
	    
	}
    }
}
