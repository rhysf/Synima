#!/usr/bin/env perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../../modules/";
#use Fasta_reader;
use read_FASTA;

# Opening statement
my $usage = "Usage: perl $0 <pep_malign.mfa> <seqs.cds>\n";
die $usage unless(@ARGV eq 2);
my ($pep_malign_file, $cds_file) = @ARGV;

# save FASTA files
my $pep_aligns = fastafile::fasta_to_struct($pep_malign_file);
my $cds_seqs = fastafile::fasta_to_struct($cds_file);

# Go through and add dashes for CDS where dashes found in peptide alignment
foreach my $acc (sort keys %{$$pep_aligns{'seq'}}) {
	my $align = $$pep_aligns{'seq'}{$acc};
	my @align_chars = split(//, $align);
	my $cds_seq = $$cds_seqs{'seq'}{$acc};
	my @cds_chars = split(//, $cds_seq);
	my $codon_align = "";
	foreach my $align_char (@align_chars) {
		if($align_char =~ /\w/) {
			my @codon_chars;
			for (1..3) {
				push (@codon_chars, shift @cds_chars);
			}
			my $codon = join("", @codon_chars);

			$codon_align .= $codon;
		}
		else { $codon_align .= "---"; }
	}
	print ">$acc\n$codon_align\n";
}
