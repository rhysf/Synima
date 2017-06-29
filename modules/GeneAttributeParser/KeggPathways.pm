package GeneAttributeParser::KeggPathways;

use strict;
use warnings;
use Carp;
use base qw (GeneAttributeParser::BaseParser);

my %kegg_pathway_definitions;

sub new {
	my $packagename = shift;

	my $self = { gene_to_KeggPathways => {},

			 };

	bless ($self, $packagename);

	return($self);

}

sub parse_dump_file {
	my $self = shift;
	my ($file) = @_;

	open (my $fh, $file) or croak ("Error, cannot open file $file");
	while (<$fh>) {
        #print;
		chomp;
		my @x = split(/\t/);

		my $gene_id = $x[1];
		my $ko_num = $x[5];
		my $definition = $x[6];
		
		if ( not defined $definition ){
			die "Error: missing KO (column 7) from file $file\n";
		}
        
		if (exists $kegg_pathway_definitions{$ko_num}) {
            
			unless ($kegg_pathway_definitions{$ko_num} eq $definition) {
				croak "$file has conflicting definitions for Ko: $ko_num GeneID: $gene_id\nNew:\'"
					. $definition . "\'\nvs\nPre-stored:\'" . $kegg_pathway_definitions{$ko_num} . "\'";
            }
        }
        
        else {
			$kegg_pathway_definitions{$ko_num} = $definition;
		}
		
		## assign EC to gene
		$self->{gene_to_KeggPathways}->{$gene_id}->{$ko_num} = 1;

        #print "$gene_id => $ko_num => $definition\n";


	}
	close $fh;
	
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
	if (my $att_href = $self->{gene_to_KeggPathways}->{$gene_id}) {

		my $annot = "";
		foreach my $ko (keys %$att_href) {
			if ($annot) {
				$annot .= ";";
			}
					
			$annot .= "$ko" . "[" . $kegg_pathway_definitions{$ko} . "]";
		}
        
		return($annot);
		
	}
	
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;

	my @genes = keys %{$self->{gene_to_KeggPathways}};

	return(@genes);
}


sub toString {
    my $self = shift;
    
    my $text = "";
    
    my @genes = $self->get_genes_with_annotations();
    foreach my $gene (@genes) {
        my $annot = $self->get_annotation($gene);
        $text .= "$gene\t$annot\n";
    }
    
    
    return($text);
}



1; #EOM

