package GeneAttributeParser::BaseParser;

## Abstract class defining required methods for superclasses

use strict;
use warnings;
use Carp;

sub new {
	
	croak "abstract class";
}
	

sub parse_dump_file {
	croak "abstract class";
}


sub get_annotation {
	croak "abstract class";
}


sub get_genes_with_annotations {
	croak "abstract class";
}



1;  #EOM

