package DataSpecFileParser;
use strict;
use warnings;
use Carp;
use FindBin;
use File::Basename;
use Data::Dumper;

#########################################
## Standard attribute types include:
#
#  Genome 
#  Annotation 
#  PFAM 
#  TIGRfam 
#  KEGG    
#  KEGG_EC 
#  COG 
#  GO_file 
#  TMHMM 
#  SignalP 
#
#########################################

our %ATT_TYPES_TO_FILE_EXTENSIONS = (	'Genome' => "genome.fa",
                                	'Annotation' => "annotation.gff3",
					'AnnotationBasic' => "annotation",
					'PFAM' => "Pfam",
                                     	'TIGRfam' => "TIGRfam",
                                     	'KEGG' => "Kegg_Blast",
                                     	'KEGG_EC' => "Kegg_EC",
                                     	'COG' => "COG",
                                     	'GO_file' => "annotation.GO",
                                     	'TMHMM' => "TMHMM",
                                     	'SignalP' => "SignalP",
                                     	'SigP4' => "SigP4",
                                     	'ARDB_file' => "ARDB",
                                     	"TBLASTN" => "alignment.gff3",
                                     	"BLAT_EST" => "alignment.gff3",
                                     	"GENEWISE" => "prediction.gff3",
                                     	"AUGUSTUS" => "prediction.gff3",
                                     	"GENEMARK" => "prediction.gff3",
                                     	"GLIMMERHMM" => "prediction.gff3",
                                     	"METAGENE" => "prediction.gff3",
                                     	"PRODIGAL" => "prediction.gff3",
                                     	"GLIMMER" => "prediction.gff3",
                                     	'RNA' => 'rna',
                                     	'TRNA' => 'trna',
                                     	'RNA_GFF' => 'rna.gff3',
                                     	'TRNA_GFF' => 'trna.gff3',
                                     
					## misc / derived
                                	'PEP' => 'annotation.pep',
                                     	'CDS' => 'annotation.cds',                                                                          
                                     	'KEGG_PATHWAYS' => 'Kegg_Blast.Kegg_Pathways',
);
my @GENE_PREDICTION_TYPES = qw(GENEWISE AUGUSTUS GENEMARK GLIMMERHMM METAGENE PRODIGAL GLIMMER);
my @GENOME_ALIGNMENT_TYPES = qw (TBLASTN BLAT_EST);

sub new {
	my $packagename = shift;
	my ($data_spec_file) = @_;

	# Is the data repo file size=0
	confess "Error, cannot find file $data_spec_file\n" unless(-s $data_spec_file);

	# Directory of datas repo file
	my $repo_dir = dirname($data_spec_file);

	# Create structure containing filename, directory, and data
	my $self = { 
		loc_data => {},
		data_spec_filename => $data_spec_file,
		repo_dir => $repo_dir,
	};

	# Make data structure an object of the package name class
	bless ($self, $packagename);
	$self->_parseDataSpecFile($data_spec_file);
	return($self);
}

sub get_genome_list {
	my $self = shift;
	my @genomes = keys %{$self->{loc_data}};
	return(@genomes);
}

sub get_repo_dir {
    my $self = shift;
    return($self->{repo_dir});
}

sub get_attribute_types {
    my $self = shift;
    my ($genome) = shift;
    unless ($genome) {
        confess "Error, require genome as parameter";
    }
    
    my $genome_data_href = $self->{loc_data}->{$genome};

    unless (ref $genome_data_href) {
        confess "Error, genome: $genome is not stored in loc data";
    }
    
    my @att_types = keys %$genome_data_href;
    
    return(@att_types);
}

sub get_attribute_value {
    my $self = shift;
    my ($genome, $att_type) = @_;

    unless (defined $genome && defined $att_type) {
        confess "Error, require params: genome, att_type";
    }

    unless (exists $self->{loc_data}->{$genome}->{$att_type}) {
        confess "Error, genome -> $att_type has no value";
    }

    return($self->{loc_data}->{$genome}->{$att_type});
}

sub get_genome_repo_dir {
	my $self = shift;
	my ($genome) = @_;
	confess "Error, need genome as param" unless($genome);
	my $genome_dir = $self->{repo_dir} . "/$genome";
	return($genome_dir);
}

sub get_data_dump_filename {
	my $self = shift;
	my ($genome, $att_type) = @_;

	my $att_type_for_filename = $att_type;

	# check for misc/derived types.
	# use the Annotation att type for base name
	if($att_type =~ m/PEP|CDS|AnnotationBasic|GO_file|SigP4|ARDB_file/) {
		$att_type_for_filename = "Annotation";
	} elsif ($att_type eq "KEGG_PATHWAYS") {
		$att_type_for_filename = "KEGG";
	} elsif( $att_type eq "RNA_GFF" ){
		$att_type_for_filename = "RNA";
	} elsif( $att_type eq "TRNA_GFF" ){
		$att_type_for_filename = "TRNA";
	}
	my $att_val_aref = $self->get_attribute_value($genome, $att_type_for_filename);

	my @filenames;
	foreach my $analysisRun (@$att_val_aref) {
		my $filename = $self->get_data_dump_filename_for_analysisRun($genome, $att_type, $analysisRun);
		#print "Filename: $filename\n";
		push (@filenames, $filename);
	}

        # Want array
	if (wantarray) { return (@filenames); } 
	# Not wanting array, just first entry
	else { return($filenames[0]); }
}

sub get_data_dump_filename_for_analysisRun {
	my $self = shift;
	my ($genome, $att_type, $analysisRun) = @_;
	confess "require analysisRun as parameter value" unless ($analysisRun);

	my $extension = $ATT_TYPES_TO_FILE_EXTENSIONS{$att_type};
	my $att_val = $analysisRun;

	if ($att_val =~ m|/|) {
		# could be path to file
		$att_val = basename($att_val);
	}

	my $filename = $self->{repo_dir} . "/$genome/$att_val" . "." . $extension;
	#print STDERR "Filename: $genome, $att_type, $analysisRun => $filename\n";
	return($filename);
}

sub get_gene_prediction_types {
	my $self = shift;
	return (@GENE_PREDICTION_TYPES);
}

sub get_genome_alignment_types {
	my $self = shift;
	return (@GENOME_ALIGNMENT_TYPES);
}

sub _parseDataSpecFile {
	my $self = shift;
	my ($data_spec_file) = @_;

	# Init loc_data
	my %genome_to_data;
	my %loc_data = ();

	# Open data repo file
	open (my $fh, $data_spec_file) or die "Error, cannot open file $data_spec_file: $!\n";
	REPO: while (my $line=<$fh>) {
		next REPO unless ($line =~ m/\S/);
		chomp $line;

		# Begining or End of section
		if(($line =~ m|^\\\\|) || ($line =~ m|^//|)) {
			if (%loc_data) {
				&_process_data(\%loc_data, \%genome_to_data);
				%loc_data = ();
			}
		} else {
			my ($key, $val) = split(/\s+/, $line);
			if ((defined $val) && ($val ne ".")) {
				push (@{$loc_data{$key}}, $val);
			}
		}
	}
	close $fh;

	# Anything left?
	if (%loc_data) { &_process_data(\%loc_data, \%genome_to_data); }

	# print STDERR Dumper(\%genome_to_data);
	$self->{loc_data} = \%genome_to_data;
	return;
}

sub _process_data {
	my ($loc_data_href, $genome_to_data_href) = @_;

	# only one genome entry allowed here.
	my $genome = $loc_data_href->{Genome}->[0]; 
	unless ($genome) { die "Error, no genome defined: " . Dumper($loc_data_href); }

	# Store everything in genome_to_data_href
	foreach my $key (keys %$loc_data_href) {
		$genome_to_data_href->{$genome}->{$key} = $loc_data_href->{$key};
	}
	return;
}

1; #EOM
