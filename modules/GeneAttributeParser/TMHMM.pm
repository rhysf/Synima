package GeneAttributeParser::TMHMM;

use strict;
use warnings;
use Carp;
use base qw (GeneAttributeParser::BaseParser);


sub new {
	my $packagename = shift;

	my $self = { gene_to_TMHMM => {},

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
		my @x = split(/\t/);

		my $gene_id = $x[1];
		my $TMHMM_count = $x[6];
		
		$self->{gene_to_TMHMM}->{$gene_id} = $TMHMM_count;
		

	}
	close $fh;
	
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
	if (my $TMHMM_count = $self->{gene_to_TMHMM}->{$gene_id}) {
		
		return("TMH:$TMHMM_count");
	}
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;

	my @genes = keys %{$self->{gene_to_TMHMM}};

	return(@genes);
}


1; #EOM

