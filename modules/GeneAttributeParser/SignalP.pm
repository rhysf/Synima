package GeneAttributeParser::SignalP;
use strict;
use warnings;
use Carp;
use base qw (GeneAttributeParser::BaseParser);

sub new {
	my $packagename = shift;
	my $self = { gene_to_SignalP => {},
			 };
	bless ($self, $packagename);
	return($self);
}

sub parse_dump_file {
	my $self = shift;
	my ($file) = @_;
	open (my $fh, $file) or croak ("Error, cannot open file $file");
	while (<$fh>) {
		chomp;
		my ($trans_id, $gene_id, $locus, $length, $predictionCode, $SProb, $SprobSignal) = split(/\t/);
		$self->{gene_to_SignalP}->{$gene_id} = { trans_id => $trans_id,
							length => $length,
							code => $predictionCode,
							SProb => $SProb,
							SProbSignal => $SprobSignal,
							};
	}
	close $fh;
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
	if (my $sigP_ref = $self->{gene_to_SignalP}->{$gene_id}) {
		my $annot = "SigP[len:" . $sigP_ref->{length} . ",SProb:" . $sigP_ref->{SProb} . "]";;
		return($annot);
	}
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;
	my @genes = keys %{$self->{gene_to_SignalP}};
	return(@genes);
}


1; #EOM

