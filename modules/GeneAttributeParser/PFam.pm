package GeneAttributeParser::PFam;
use strict;
use warnings;
use Carp qw (confess cluck croak);
use base qw (GeneAttributeParser::BaseParser);

my %PFAM_definitions;

sub new {
	my $packagename = shift;
	my ($cutoff) = @_;  ## must be 'TRUSTED' or 'NOISE'
	
	unless ($cutoff eq "TRUSTED" || $cutoff eq "NOISE") {
		croak "must specify score cutoff as TRUSTED or NOISE";
	}
	
	my $self = { gene_to_PFAM => {},
				 cutoff => $cutoff,
			 };

	bless ($self, $packagename);

	return($self);

}

sub parse_dump_file {
	my $self = shift;
	my ($file) = @_;

	open (my $fh, $file) or croak ("Error, cannot open file $file");
	while (<$fh>) {
		my $line = $_;
		if (/^\#/) { next; }
		chomp;
		my @x = split(/\t/);
		
		my $gene_id = $x[1];
		my $pfam_acc = $x[9];
		my $pfam_score = $x[15];
		my $hmm_acc = $x[16];
		my $pfam_descr = $x[18];
		my $noise_cutoff = $x[19];
		my $trusted_cutoff = $x[20];
        

        unless ($pfam_score =~ /\d/ && $trusted_cutoff =~ /\d/) {
            print STDERR "-warning, $file contains line: $line with uninterpretable scores; ignoring it.\n";
            next;
        }
        
        
		if (  
			  ($self->{cutoff} eq "NOISE" && $noise_cutoff =~ /\d/ && $pfam_score >= $noise_cutoff)
			  ||
			  ($self->{cutoff} eq "TRUSTED" && $pfam_score >= $trusted_cutoff) 
			  )
		{
			
			
			unless ($hmm_acc && $pfam_acc && $pfam_descr) {
				print STDERR "Line of file: $file: $line lacks hmm accession or pfam accession or description.\n";
				next;
			}
			
			if ($hmm_acc ne $pfam_acc) {
				$pfam_acc = "$hmm_acc|$pfam_acc"; # bundle
			}
			
			## assign pfam to gene
			$self->{gene_to_PFAM}->{$gene_id}->{$pfam_acc} = 1;
			
			if (exists $PFAM_definitions{$pfam_acc} &&
				$PFAM_definitions{$pfam_acc} ne $pfam_descr) {
			    cluck "Error, pfam accession: $pfam_acc with descr [$pfam_descr]\n"
					. "is not equivalent to previously stored value of:\n"
					. "[$PFAM_definitions{$pfam_acc}]  ";
			}
			else {
				# store definition
				$PFAM_definitions{$pfam_acc} = $pfam_descr;
			}
		}
	}
	close $fh;
	return;
}

sub get_annotation {
	my $self = shift;
	my ($gene_id) = @_;
	
	if (my $att_href = $self->{gene_to_PFAM}->{$gene_id}) {
		
		my $annot = "";
		foreach my $pfam_acc (keys %$att_href) {
			if ($annot) {
				$annot .= ";";
			}
			my $descr = $PFAM_definitions{$pfam_acc} || "none";
			$annot .= "$pfam_acc" . "[" . $descr . "]";
		}
		return($annot);
	}
	else {
		return("");
	}
}

sub get_genes_with_annotations {
	my $self = shift;
	my @genes = keys (%{$self->{gene_to_PFAM}});
	return(@genes);
}

1; #EOM
