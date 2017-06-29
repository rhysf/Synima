package GeneAttributeParser::GO;
use strict;
use warnings;
use Carp;
use base qw (GeneAttributeParser::BaseParser);

sub new {
	my $packagename = shift;
	my $self = { gene_to_GO => {},
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
		my $gene_info_txt = $x[0];
		my $GO_id = $x[1];
		
        my (@gene_id_line) = split(/\s+/, $gene_info_txt);
        my $gene_id = $gene_id_line[1];
        $gene_id =~ s/gene_id=//g;

		$self->{gene_to_GO}->{$gene_id}->{$GO_id} = 1;
        #print STDERR "$gene_id => $GO_id\n";
        
	}
	close $fh;
	
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
    unless ($gene_id) {
        confess "Error, need gene_id as parameter";
    }
    
	if (my $att_href = $self->{gene_to_GO}->{$gene_id}) {
		
		my $annot = "";
		foreach my $GO_id (keys %$att_href) {
			if ($annot) {
				$annot .= ";";
			}
			$annot .= "$GO_id";
		}
		return($annot);
	}
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;
	my @genes = keys %{$self->{gene_to_GO}};
	return(@genes);
}

sub toString {
    my $self = shift;
    my $text = "";
    my @genes = $self->get_genes_with_annotations();
    
    foreach my $gene (@genes) {
        $text .= "$gene\t" . $self->get_annotation($gene) . "\n";
    }
    return($text);
}

1; #EOM
