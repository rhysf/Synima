#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use Synima;
use read_GFF;
use read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <Repo_spec.txt>
Optional: -f Feature wanted from GFF [mRNA]
	  -s Seperator in GFF description for gene names (\" ; etc) [;]
	  -d GFF description part number with the parent/gene info [0]
	  -m Remove additional comments in column [ID=]
Notes: Will copy all transcripts and specified features from GFF into primary fasta files\n";
our($opt_r, $opt_f, $opt_s, $opt_d, $opt_m);
getopt('rfsdm');
die $usage unless ($opt_r);
die "Cannot open $opt_r : $!\n" unless (-e $opt_r);
if(!defined $opt_f) { $opt_f = 'mRNA'; }
if(!defined $opt_s) { $opt_s = ';'; }
if(!defined $opt_d) { $opt_d = 0; }
if(!defined $opt_m) { $opt_m = 'ID='; }

# Perform data retrievals
my $data_manager = new DataSpecFileParser($opt_r);
my @genomes = $data_manager->get_genome_list();
foreach my $genome (@genomes) {
	warn "Indexing $genome...\n";		
	my $genome_fasta = $data_manager->get_data_dump_filename($genome, 'Genome');
}

# Create repository sequence databases
&create_full_repo_sequence_databases($data_manager, \@genomes, $opt_r);

# Get all annotation files together into one gff3
my $all_annotations = "$opt_r.all.GFF3";
gfffile::combine_all_gff3_files_in_repo($opt_r, $all_annotations, $opt_f, $opt_s, $opt_d, $opt_m);

# Check FASTA and GFFs match
warn "Checking FASTA and GFF repository sequence databases match...\n";
foreach my $type(qw(PEP CDS)) {
	my $fasta = "$opt_r.all.$type";
	my ($sequences, $descriptions, $order) = fastafile::fasta_id_to_seq_hash($fasta);
	my ($gff, $strand) = gfffile::gff_to_contig_parent_to_cds_hash($all_annotations);
	my ($found, $not_found, $total) = (0, 0, 0);
	foreach my $contig(keys %{$gff}) {
		foreach my $parent(keys %{$$gff{$contig}}) {
			if(defined $$sequences{$parent}) { $found++; }
			else { $not_found++; }
		}
	}
	$total = ($found + $not_found);
	warn "\n\n\n$found / $total found\n";
	if($found eq $total) { warn "$fasta and $all_annotations repository sequence databases are correctly formatted.\n"; }
	else { warn "$fasta and $all_annotations repository sequence databases are not correctly formatted. Change settings and re-run, or rename ID's in FASTA or GFF to match.\n"; }
}
warn "Finished.\n";

sub create_full_repo_sequence_databases {
	my ($data_manager, $genomes_aref, $repo) = @_;

	# Copy across CDS and PEP for each genome
	warn "Creating repository sequence databases...\n";
	foreach my $genome (@$genomes_aref) {
		foreach my $type(qw(PEP CDS)) {
			my $file = $data_manager->get_data_dump_filename($genome, $type) or die "Error, cannot find $type file for genome: $genome : $!\n";
			my $ofile = "$repo.all.$type";
			warn "Copying $file to $ofile...\n";
			&parse_and_join_fasta_to_outfiles($file, $ofile, $genome);
		}
	}
	return 1;
}

sub parse_and_join_fasta_to_outfiles {
	my ($file, $ofile, $genome) = @_;
	open my $fh,  '<', $file  or die "Error, cannot open file $file : $!\n";
	open my $ofh, '>>', $ofile or die "Error, cannot write to $ofile : $!\n";
	while (my $line=<$fh>) { 
		chomp $line;
		$line = &make_broad_style_id_for_fasta($line, $genome);
		print $ofh "$line\n"; 
	}
	close $fh;
	close $ofh;
	return 1;
}

# Broad Format = >7000011728610201 gene_id=7000011728610200 locus=None name="flagellum-specific ATP synthase" genome=Esch_coli_MGH121_V1 analysisRun=Esch_coli_MGH121_V1_POSTPRODIGAL_2	
sub make_broad_style_id_for_fasta {
	my ($line, $genome) = @_;
	if(($line =~ m/^>/) && (($line !~ m/gene_id\=/) && ($line !~ m/locus\=/) && ($line !~ m/name\=/) && ($line !~ m/genome\=/) && ($line !~ m/analysisRun\=/))) { 
		my @id_parts = split / /, $line;
		my $id = $id_parts[0];
		$id =~ s/^>//;
		$line = ">$id gene_id=$id locus=$id name=\"hypothetical protein\" genome=$genome analysisRun=$genome\_1";
	}
	return $line;
}
