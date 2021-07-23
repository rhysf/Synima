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

	# Remove carriage returns
	$line =~ s/\r//g;

	# Check ID line for all Broad-style parts. Adding them if not present
	if($line =~ m/^>/) {
		my $id;
		my $new_line;
		#warn "old id = $line\n";
		my @id_parts = split / /, $line;
		for(my $i=0; $i<6; $i++) {
			my $id_part = '';
			if(defined $id_parts[$i]) { $id_part = $id_parts[$i]; }

			# ID
			if($i eq 0) { 
				my @id_parts2 = split />/, $id_part;
				$id = $id_parts[1];
				$new_line = $id_part; 
			}
			# Gene id
			elsif($i eq 1) {
				if($id_part =~ m/gene_id\=/) { $new_line .= " $id_part"; }
				else { $new_line .= " gene_id=$id"; }
			}
			# Locus
			elsif($i eq 2) {
				if($id_part =~ m/locus\=/) { $new_line .= " $id_part"; }
				else { $new_line .= " locus=$id"; }
			}
			# Name
			elsif($i eq 3) {
				if($id_part =~ m/name\=\"/) { 
					# spaces in name
					my $name;
					my $remove_extra_name_elements_count = 0;
					NAME_PARTS: for(my $j=3; $j <= scalar(@id_parts); $j++) {
						$name .= " $id_parts[$j]";
						#warn "j $j = $name\n";
						# end of name?
						if($id_parts[$j] =~ m/\"$/) {
							#warn "j $j end of name\n";
							if($j > 3) { $remove_extra_name_elements_count++; }
							last NAME_PARTS;
						}

						# Remove extra id parts from id array (name should fill only 1 element)
						if($j > 3) { $remove_extra_name_elements_count++; }
					}

					# remove name parts from rest of id array (name should fill only 1 element)
					REMOVE_NAME_PARTS: for(my $k=0; $k<$remove_extra_name_elements_count; $k++) {
						#warn "Remove $id_parts[4]\n";
						splice @id_parts, 4, 1; 
					}
					$new_line .= " $name"; 
				}
				else { $new_line .= " name=\"hypothetical protein\""; }
			}
			# Genome
			elsif($i eq 4) {
				if($id_part =~ m/genome\=/) { $new_line .= " $id_part"; }
				else { $new_line .= " genome=$genome"; }
			}
			# AnalysisRun
			elsif($i eq 5) {
				if($id_part =~ m/analysisRun\=/) { $new_line .= " $id_part"; }
				else { $new_line .= " analysisRun=$genome\_1"; }
			}
			#warn "$i) $id_part\n";
		}
		#warn "new id = $new_line\n\n";
		$line = $new_line;
	}

	return $line;
}
