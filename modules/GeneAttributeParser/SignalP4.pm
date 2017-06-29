package GeneAttributeParser::SignalP4;

use strict;
use warnings;
use Carp;
use base qw (GeneAttributeParser::BaseParser);


sub new {
	my $packagename = shift;

	my $self = { gene_to_SignalP4 => {},

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

        ##need to ignore comment lines
        
        my @x = split(/\s+/);
        if (exists $x[9]) {
            if ($x[9] eq 'Y') {
                my $trans_id = $x[0];
                my $cleave_site = $x[4];
                my $D_score = $x[8];

                my $gene_id = $x[0]; ###### $x[0] is actually transcript id, gene id not in output
		
		        $self->{gene_to_SignalP4}->{$gene_id} = { trans_id => $trans_id,
												         cleave => $cleave_site,
												         dscore => $D_score,
											     };
			}
		
	    }
	}
	close $fh;
	
	
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
	if (my $sigP4_ref = $self->{gene_to_SignalP4}->{$gene_id}) {
		
		my $annot = "SigP4[DScore:" . $sigP4_ref->{dscore} . "]";;
		
		return($annot);
	}
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;
   
	my @genes = keys %{$self->{gene_to_SignalP4}};

	return(@genes);
}


1; #EOM

