package GeneAttributeParser::KeggBlast;

use strict;
use warnings;
use Carp;
use base qw (GeneAttributeParser::BaseParser);

my %kegg_top_hits;
my %kegg_top_scores;

sub new {
	my $packagename = shift;

	my $self = { gene_to_KeggBlast => {},

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
		if (/^#/) {
		} else {
		    my @x = split(/\t/);

		    my $gene_id = $x[1];
		    my $hit_name = $x[4];
		    my $hit_score = $x[9];
		    
		    if ($x[4] =~ /\w+:\w+  ([a-z]{3}[A-Z|0-9]{0,2});/) {
        
		        if (exists $kegg_top_hits{$gene_id} and exists $kegg_top_scores{$gene_id}) {
                    if ($hit_score > $kegg_top_scores{$gene_id}) {
			            $kegg_top_hits{$gene_id} = $hit_name;
			            $kegg_top_scores{$gene_id} = $hit_score; 
			    
			            $self->{gene_to_KeggBlast}->{$gene_id}->{$hit_name} = 1;
                    }
                }
        
                else {
			        $kegg_top_hits{$gene_id} = $hit_name;
			        $kegg_top_scores{$gene_id} = $hit_score;
			
			        $self->{gene_to_KeggBlast}->{$gene_id}->{$hit_name} = 1;
			    }
		    }
		
		    ## assign EC to gene
		    #$self->{gene_to_KeggBlast}->{$gene_id}->{$ko_num} = 1;

            #print "$gene_id => $ko_num => $definition\n";

        }
	}
	close $fh;
	
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
	if (my $att_href = $self->{gene_to_KeggBlast}->{$gene_id}) {

		my $annot = "$kegg_top_hits{$gene_id}";

		return($annot);
		
	}
	
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;

	my @genes = keys %{$self->{gene_to_KeggBlast}};

	return(@genes);
}


#sub toString {
#    my $self = shift;
#    
#    my $text = "";
#    
#    my @genes = $self->get_genes_with_annotations();
#    foreach my $gene (@genes) {
#        my $annot = $self->get_annotation($gene);
#        $text .= "$gene\t$annot\n";
#    }
#    
#    
#    return($text);
#}



1; #EOM

