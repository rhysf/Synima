package GeneAttributeParser::KeggEC;

use strict;
use warnings;
use Carp;
use base qw (GeneAttributeParser::BaseParser);

my %EC_definitions;

sub new {
	my $packagename = shift;

	my $self = { gene_to_EC => {},

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

		my $gene_id = $x[0];
		my $EC = $x[6];
		my $definition = $x[8];
		
		$definition =~ s/^\s+//;
		$definition =~ s/\s+$//;

		if (exists $EC_definitions{$EC}) {

			unless ($EC_definitions{$EC} eq $definition) {
				croak "$file has conflicting definitions for EC: $EC GeneID: $gene_id\n New:\'"
					. $definition . "\'\nvs\nPre-stored\'" . $EC_definitions{$EC} . "\'";
			}
		}
		else {
			$EC_definitions{$EC} = $definition;
		}
		
		## assign EC to gene
		$self->{gene_to_EC}->{$gene_id}->{$EC} = 1;


	}
	close $fh;
	
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
	if (my $att_href = $self->{gene_to_EC}->{$gene_id}) {

		my $annot = "";
		foreach my $EC (keys %$att_href) {
			if ($annot) {
				$annot .= ";";
			}
					
			$annot .= "$EC" . "[" . $EC_definitions{$EC} . "]";
		}

		return($annot);
		
	}
	
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;

	my @genes = keys %{$self->{gene_to_EC}};

	return(@genes);
}


1; #EOM

