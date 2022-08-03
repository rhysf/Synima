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
          -v Verbose mode (report IDs missing) (y/n) [n]
Notes: Will QC and copy all transcripts and specified features from GFF into primary GFF and fasta files\n";
our($opt_d, $opt_f, $opt_r, $opt_s, $opt_m, $opt_v);
getopt('dfrsmv');
die $usage unless ($opt_r);
die "Cannot open $opt_r : $!\n" unless (-e $opt_r);
if(!defined $opt_d) { $opt_d = 0; }
if(!defined $opt_f) { $opt_f = 'mRNA'; }
if(!defined $opt_s) { $opt_s = ';'; }
if(!defined $opt_m) { $opt_m = 'ID='; }
if(!defined $opt_v) { $opt_v = 'n'; }

# Combine all annotation files together into one gff3
my $all_annotations = "$opt_r.all.GFF3";
# new_transcript_ids hash = genome -> old_to_new/new_to_old -> id1 = id2
my $new_transcript_ids = gfffile::combine_all_gff3_files_in_repo($opt_r, $all_annotations, $opt_f, $opt_s, $opt_d, $opt_m);

# Perform data retrievals from repo_spec
my $data_manager = new DataSpecFileParser($opt_r);
my @genomes = $data_manager->get_genome_list();
foreach my $genome (@genomes) {
	warn "Indexing $genome...\n";		
	my $genome_fasta = $data_manager->get_data_dump_filename($genome, 'Genome');
}

# Combine all CDS or PEP FASTA files for each genome
warn "Creating repository sequence databases...\n";
foreach my $genome (@genomes) {
	foreach my $type(qw(PEP CDS)) {
		my $file = $data_manager->get_data_dump_filename($genome, $type) or die "Error, cannot find $type file for genome: $genome : $!\n";
		my $ofile = "$opt_r.all.$type";
		warn "Copying $file to $ofile for $genome...\n";
		&parse_and_join_fasta_to_outfiles($file, $ofile, $genome, $new_transcript_ids);
	}
}

# Check FASTA and GFFs match
warn "Checking FASTA and GFF repository sequence databases match...\n";
foreach my $type(qw(PEP CDS)) {
	my $fasta = "$opt_r.all.$type";
	my ($sequences) = fastafile::fasta_to_struct($fasta);
	my ($gff, $strand) = gfffile::gff_to_contig_parent_to_cds_hash($all_annotations);

	# Looping through every GFF entry and looking for corresponding FASTA entries
	my %genes_in_GFF;
	my ($found, $not_found, $total) = (0, 0, 0);
	foreach my $contig(sort keys %{$gff}) {
		foreach my $parent(sort keys %{$$gff{$contig}}) {

			# Check if gene from Repo_spec.txt.all.GFF3 is in Repo_spec.txt.all.CDS/PEP
			$genes_in_GFF{$parent} = 1;
			if(defined $$sequences{'seq'}{$parent}) { $found++; }
			else { 
				$not_found++; 
				if($opt_v ne 'n') { warn "WARNING: No gene $parent from GFF found in FASTA\n"; }
			}
		}
	}
	$total = ($found + $not_found);

	# Looping through every FASTA entry and looking for corresponding GFF entries
	my ($found2, $not_found2, $total2) = (0, 0, 0);
	foreach my $gene(keys %{$$sequences{'seq'}}) {
		if(defined $genes_in_GFF{$gene}) { $found2++; }
		else {
			$not_found2++;
			if($opt_v ne 'n') { warn "WARNING: No gene $gene from FASTA found in GFF\n"; }
		}

	}
	$total2 = ($found2 + $not_found2);

	warn "\n$found / $total found (GFF $all_annotations entries in FASTA $fasta)\n";
	warn "$found2 / $total2 found (FASTA $fasta entries in GFF $all_annotations)\n";
	if(($found eq $total) && ($found2 eq $total2)) { warn "$fasta and $all_annotations repository sequence databases are correctly formatted.\n"; }
	else { warn "WARNING: $fasta and $all_annotations repository sequence databases are not correctly formatted. Change settings and re-run, or rename ID's in FASTA or GFF to match. Use -v for further info.\n"; }
}
warn "$0: finished check.\n";

sub parse_and_join_fasta_to_outfiles {
	my ($file, $ofile, $genome, $new_transcript_ids) = @_;
	open my $fh,  '<', $file  or die "Error, cannot open file $file : $!\n";
	open my $ofh, '>>', $ofile or die "Error, cannot write to $ofile : $!\n";
	while (my $line=<$fh>) { 
		chomp $line;
		$line = &make_broad_style_id_for_fasta($line, $genome, $file, $new_transcript_ids);
		print $ofh "$line\n"; 
	}
	close $fh;
	close $ofh;
	return 1;
}

# Broad Format = >7000011728610201 gene_id=7000011728610200 locus=None name="flagellum-specific ATP synthase" genome=Esch_coli_MGH121_V1 analysisRun=Esch_coli_MGH121_V1_POSTPRODIGAL_2	
sub make_broad_style_id_for_fasta {
	my ($line, $genome, $file_name, $new_transcript_ids) = @_;

	# Remove carriage returns
	$line =~ s/\r//g;
	return $line if($line !~ m/^>/);

	# AnalysisRun
	my $analysisRun = fileparse($file_name);
	$analysisRun =~ s/\.annotation\.(pep|cds)//i;

	# Check ID line for all Broad-style parts. Adding them if not present
	my $new_line;
	#warn "old id = $line\n";
	my @id_parts = split / /, $line;
	die "make_broad_style_id_for_fasta: Error: Unrecognised ID given for $line in FASTA (for $genome / $file_name)\n" if(scalar(@id_parts) < 1);

	# Transcript ID
	my @id_parts2 = split />/, $id_parts[0];
	my $id = $id_parts2[1];

	# Unique Transcript ID may have been made for GFF
	if(defined $$new_transcript_ids{$genome}{'old_to_new'}{$id}) {
		$id = $$new_transcript_ids{$genome}{'old_to_new'}{$id};
	}

	# Gene ID
	my $gene_id = "gene_id=$id";
	if(defined $id_parts[1]) {
		if($id_parts[1] =~ m/gene_id\=/) { $gene_id = $id_parts[1]; }
	}

	# Locus ID
	my $locus_id = "locus=$id";
	if(defined $id_parts[2]) {
		if($id_parts[2] =~ m/locus\=/) { $locus_id = $id_parts[2]; }
	}

	# Name
	my $name = "name=\"hypothetical protein\"";
	if(defined $id_parts[3]) {
		if($id_parts[3] =~ m/name\=\"/) { 

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
			$name = "$name"; 
		}
	}

	# Genome
	my $genome_id = "genome=$genome";
	my $analysis_run = "analysisRun=$analysisRun";

	# return new line
	$line = ">$id $gene_id $locus_id $name $genome_id $analysis_run";
	return $line;
}
