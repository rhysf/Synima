#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use Synima;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <Repo_spec.txt>
Notes: Will copy all transcripts for each gene into primary fasta files\n";
our($opt_r);
getopt('r');
die $usage unless ($opt_r);
die "Cannot open $opt_r : $!\n" unless (-e $opt_r);

# Perform data retrievals
my $data_manager = new DataSpecFileParser($opt_r);
my @genomes = $data_manager->get_genome_list();
foreach my $genome (@genomes) {
	print "// Indexing $genome\n";		
	my $genome_fasta = $data_manager->get_data_dump_filename($genome, 'Genome');
}

# Generate individual genome attribute summaries.
#print STDERR "-Summarizing attributes for each genome\n";
#my $cmd = "$misc_util/describe_gene_attributes.pl $opt_r";
#synima::process_cmd($cmd);

# Create database files for pep and cds and make blastable
my $pep_database = "$opt_r.all.pep";
my $cds_database = "$opt_r.all.cds";
&create_full_repo_sequence_databases($data_manager, \@genomes, $pep_database, $cds_database);

sub create_full_repo_sequence_databases {
	my ($data_manager, $genomes_aref, $pep_db_filename, $cds_db_filename) = @_;

	# Open output files
	open my $cds_fh, '>', $cds_db_filename or die "Error, cannot write to $cds_db_filename : $!\n";
	open my $pep_fh, '>', $pep_db_filename or die "Error, cannot write to $pep_db_filename : $!\n";

	# Copy across CDS and PEP for each genome
	# Broad Format = >7000011728610201 gene_id=7000011728610200 locus=None name="flagellum-specific ATP synthase" genome=Esch_coli_MGH121_V1 analysisRun=Esch_coli_MGH121_V1_POSTPRODIGAL_2	
	warn "Writing $cds_db_filename and $pep_db_filename...\n";
	foreach my $genome (@$genomes_aref) {
		warn "\tCopying over pep and cds for $genome\n";
		my $pep_file = $data_manager->get_data_dump_filename($genome, 'PEP') or die "Error, cannot find pep file for genome: $genome : $!\n";
		open my $fh, '<', $pep_file or die "Error, cannot open file $pep_file : $!\n";
		while (my $line=<$fh>) { 
			chomp $line;
			$line = &make_broad_style_id_for_fasta($line, $genome);
			print $pep_fh "$line\n"; 
		}
		close $fh;

	   	my $cds_file = $data_manager->get_data_dump_filename($genome, 'CDS') or die "Error, cannot find cds file for genome: $genome : $!\n";
	   	open $fh, '<', $cds_file or die "Error, cannot open $cds_file : $!\n";
	   	while (my $line=<$fh>) { 
			chomp $line;
			$line = &make_broad_style_id_for_fasta($line, $genome);
			print $cds_fh "$line\n"; 
		}
		close $fh;
	}
	warn "Finished writing $cds_db_filename and $pep_db_filename\n";

	# Close output files and return
	close $cds_fh;
	close $pep_fh;
	return 1;
}

sub make_broad_style_id_for_fasta {
	my ($line, $genome) = @_;
	if(($line =~ m/^>/) && (($line !~ m/gene_id\=/) && ($line !~ m/locus\=/) && ($line !~ m/name\=/))) { 
		my $id = $line;
		$id =~ s/^>//;

		# If there is a description - move into a single long id for parsing later
		$line =~ s/ /_/g;
		$line .= " gene_id=$id locus=$id name=\"hypothetical protein\" genome=$genome analysisRun=$genome\_1";
	}
	return $line;
}
