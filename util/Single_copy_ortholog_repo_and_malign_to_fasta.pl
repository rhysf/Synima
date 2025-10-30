#!/usr/bin/env perl -w
use strict;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use File::Basename;
use Data::Dumper;
use Getopt::Std;
use read_FASTA;

### r.farrer@exeter.ac.uk

# Opening statement
my $usage =  "Usage: perl $0 -r <Repo_spec.txt> -g <orthologAnalysisRun.clusters (ortholog data dump output)>
Optional: -c\tUse codon alignments instead of protein alignments (y/n) []
          -a\tColumn from annotation files to use for ids [0]
          -e\tRestrict to genomes (file containing the list of genomes to restrict alignments and trees to.) []
	  -t\tMake trees (y/n) [y]\n";
our($opt_a, $opt_c, $opt_e, $opt_g, $opt_r, $opt_t);
getopt('acegrt');
die $usage unless ($opt_r && $opt_g);
if(!defined $opt_a) { $opt_a = 0; }
if(!defined $opt_t) { $opt_t = 'y'; }

# Programs
my $fastTree_prog = "FastTree";
my $MFA_block_extractor = "$FindBin::Bin/support_scripts/MFA_block_extractor.pl";
unless (-e $MFA_block_extractor) { die "Error, cannot find support script $MFA_block_extractor"; }

# Find the multiple alignment folder
my $malign_dir = basename($opt_g) . ".MALIGN_DIR";
if ($opt_e) { $malign_dir = basename($opt_e) . ".MALIGN_DIR"; }
unless (-d $malign_dir) { die "Error, cannot find directory $malign_dir"; }

# Save Repo spec
my $data_manager = new DataSpecFileParser($opt_r);

# Save genome IDs we're interested in
my $genomes = &data_manager_save_genomes_of_interest($data_manager, $opt_e);

# Save transcript -> genome hash from peptide files
my %trans_id_to_genome;
foreach my $genome(@{$genomes}) {
	warn "Saving transcripts -> genome ($genome)\n";
	my $protein_file_base = $data_manager->get_data_dump_filename($genome, 'PEP');
	my $protein_file = "$protein_file_base.synima-parsed.PEP";
	my $protein_seqs = fastafile::fasta_to_struct($protein_file);
	foreach my $id(sort keys %{$$protein_seqs{'seq'}}) {
		$trans_id_to_genome{$id} = $genome;
	}
	
}

# Get list of mfa files
my $core_sc_gene_clusters = &get_single_copy_core($opt_g, $genomes);
my $alignment_files = &save_muscle_fastas_from_folder($core_sc_gene_clusters, $malign_dir);

# Concatenate alignments
my $concat_alignments_by_genome = &concatenate_alignments_by_genomes($alignment_files, \%trans_id_to_genome);

# outfile name
my $concat_align_file = "SC_core_concat.aa-based.$$.mfa";
if($opt_c) { $concat_align_file = "SC_core_concat.codon-based.$$.mfa"; }

# print
warn "Contactenating alignments into $concat_align_file...\n";
open (my $ofh, ">$concat_align_file") or die "Error, cannot write file $concat_align_file";
foreach my $acc (sort keys %{$concat_alignments_by_genome}) {
	print $ofh ">$acc\n";
	my $alignment = $$concat_alignments_by_genome{$acc};
	$alignment =~ s/(\S{60})/$1\n/g;
	chomp $alignment;
	print $ofh $alignment . "\n";
}
close $ofh;

# extract alignment blocks:
my $block_extracted_mfa_file = ($opt_c) ? "SC_core_concat.blocks_extracted.codon-based.$$.mfa" : "SC_core_concat.blocks_extracted.aa-based.$$.mfa";
my $cmd = "perl $MFA_block_extractor $concat_align_file > $block_extracted_mfa_file";
&process_cmd($cmd);


#if ($opt_c) { $fastTree_prog .= " -nt "; }

# Generate FastTree trees
#if($opt_t eq 'y') {
#	eval {
#		$cmd = "$fastTree_prog $concat_align_file > $concat_align_file.tree"; # sometimes this step fails if it's too large.
#		&process_cmd($cmd);
#	};
#	if ($@) { print "Error:  $@"; }

#	$cmd = "$fastTree_prog $block_extracted_mfa_file > $block_extracted_mfa_file.tree";
#	&process_cmd($cmd);
#}

sub process_cmd {
	my ($cmd) = @_;
	warn "CMD: $cmd\n";
	my $ret = system($cmd);
	if ($ret) { die "Error, cmd: $cmd died with ret: $ret"; }
	return;
}

sub get_single_copy_core {
	my ($gene_clusters_file, $genomes_aref) = @_;

	# Move genome ID array to hash
	my %genomes = map { + $_ => 1 } @$genomes_aref;
	my $num_genomes = scalar(@$genomes_aref);

	my %core_cluster_info;
	my $curr_cluster_id = "";
	my %cluster_to_data;

	open (my $fh, '<', $gene_clusters_file) or die "Error, cannot open file $gene_clusters_file";
	CLUSTERS: while (my $line=<$fh>) {
		next CLUSTERS unless ($line =~ m/\w/);
		my @bits = split /\t/, $line;
		my ($cluster_id, $genome, $gene_set_name, $trans_id, $gene_id, $locus, @rest) = @bits;
		next CLUSTERS unless ($genomes{$genome});

		if ($cluster_id ne $curr_cluster_id) {
			if (%cluster_to_data) {
				&process_cluster_info($curr_cluster_id, \%cluster_to_data, $num_genomes, \%core_cluster_info);
			}
			# reinit
			%cluster_to_data = ();

		}

		push (@{$cluster_to_data{$genome}}, $trans_id);
		$curr_cluster_id = $cluster_id;
	}

	# get last one
	if (%cluster_to_data) {
		&process_cluster_info($curr_cluster_id, \%cluster_to_data, $num_genomes, \%core_cluster_info);
	}
	return \%core_cluster_info;
}

sub save_muscle_fastas_from_folder {
	warn "save_muscle_fastas_from_folder: $malign_dir\n";
	my ($core_sc_gene_clusters, $malign_dir) = @_;
	my @alignment_files;
	foreach my $cluster (keys %{$core_sc_gene_clusters}) {
		my $alignment_file = "$malign_dir/$cluster.pep.mfa";
		if($opt_c) { $alignment_file = "$malign_dir/$cluster.codon.mfa"; } 
		push (@alignment_files, $alignment_file);
	}
	print "Single Copy Core consists of " . scalar(@alignment_files) . " genes.\n";
	die "Found no orthologs\n" if(scalar(@alignment_files) eq 0);
	return \@alignment_files;
}


sub process_cluster_info {
	my ($curr_cluster_id, $cluster_to_data_href, $num_genomes, $core_cluster_info_href) = @_;

	unless (scalar (keys %$cluster_to_data_href) eq $num_genomes) {
		# not core
		return;
	}
	
	# check for single copy core
	foreach my $genome (keys %$cluster_to_data_href) {
		if (scalar @{$cluster_to_data_href->{$genome}} != 1) {
			return;
			# too many genes
		}
	}

	## if made it here, got a single copy core.  Copy info over
	foreach my $genome (keys %$cluster_to_data_href) {
		
		my $trans_id = $cluster_to_data_href->{$genome}->[0];
		
		$core_cluster_info_href->{$curr_cluster_id}->{$genome} = $trans_id;
	}
	return;
}

sub concatenate_alignments_by_genomes {
	my ($alignment_files_aref, $trans_id_to_genome_href) = @_;

	warn "concatenate_alignments_by_genomes...\n";

	# Concatenate all alignments (genome -> concatenated sequence)
	my %concat_alignment;
	foreach my $alignment_file (@$alignment_files_aref) {
		
		#my $fasta_reader = new Fasta_reader($alignment_file);
		#my %seqs = $fasta_reader->retrieve_all_seqs_hash();

		my $seqs = fastafile::fasta_to_struct($alignment_file);

		warn "concatenate_alignments_by_genomes: $alignment_file\n";
		my $length_of_seq;
		foreach my $acc (sort keys %{$$seqs{'seq'}}) {
			my $genome = $$trans_id_to_genome_href{$acc} or die "Error, cannot identify genome corresponding to acc [$acc]";
			
			my $alignment = $$seqs{'seq'}{$acc};
			$concat_alignment{$genome} .= $alignment;
			#warn "genome $genome -> $acc =\n$alignment\n";
		
			# QC (different lengths of individual sequences)
			if(!defined $length_of_seq) { $length_of_seq = length($alignment); }
			else { 
				if($length_of_seq ne length($alignment)) {
					die "error: $alignment_file $genome -> $acc has different length\n";
				}
			}

		}
		#die "end here\n";
	
		# QC check concat lengths
		my $length_found;
		foreach my $genome(sort keys %concat_alignment) {
			my $seq = $concat_alignment{$genome};
			my $length = length($seq);
			if(!defined $length_found) { $length_found = $length; }
			else {
				die "error: $alignment_file $genome -> has different concat lengths: $length ne $length_found\n" if($length ne $length_found);
			}
		}
	}

	# Ensure consistent lengths
	my @align_lengths;
	my $len = undef;
	my $len_inconsistent_flag = 0;
	foreach my $acc (keys %concat_alignment) {
		my $alignment = $concat_alignment{$acc};
		my $length = length($alignment);
		push (@align_lengths, [$acc, $length]);
		
		if ($len) {
			if ($length != $len) {
				$len_inconsistent_flag = 1;
			}
		}
		else { $len = $length; }
	}
	if ($len_inconsistent_flag) {
		die "Error, lengths inconsistent: " . Dumper(\@align_lengths);
	}
	return \%concat_alignment;
}

sub data_manager_save_genomes_of_interest {
	my ($data_manager, $restrict_optional) = @_;
	my @genomes;

	# restrict
	if ($restrict_optional) {
		@genomes = `cat $restrict_optional`;
		foreach my $genome (@genomes) { $genome =~ s/\s+//g; }

		# Ignore empty entries
		@genomes = grep { /\w/ } @genomes; 
	} else { @genomes = $data_manager->get_genome_list(); }
	#warn Dumper(@genomes);
	return \@genomes;
}

