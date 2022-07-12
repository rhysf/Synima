#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use File::Basename;
use Synima;
use read_GFF;

### rfarrer@broadinstitute.org

# Opening statement
my $usage = "Usage: perl $0 -r <repo_spec> -c <Ortholog cluster data (E.g. ORTHOMCLBLASTFILE.clusters)>\n
Optional: -z File containing a list of genomes to restrict the analysis to []
	  -i Minimum number of paired genes in a single dagchain [4]
	  -o Cmds outdir [dagchainer_rundir]
	  -l Cmds outfile [cluster_cmds]
	  -g Run commands on the grid (y/n) [n]
	  -p Platform (UGER, LSF, GridEngine) [UGER]
	  -q Queue (hour, short, long) [short]
	  -v Verbose (y/n) [n]\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_e, $opt_f, $opt_g, $opt_h, $opt_i, $opt_j, $opt_k, $opt_l, $opt_m, $opt_n, $opt_o, $opt_p, $opt_q, $opt_r, $opt_s, $opt_t, $opt_u, $opt_v, $opt_w, $opt_x, $opt_y, $opt_z);
getopt('abcdefghijklmnopqrstuvwxyz');
die $usage unless ($opt_r && $opt_c);
foreach($opt_r, $opt_c) { die "Cannot find $_\n" unless(-e $_); }
if(!defined $opt_z) { $opt_z = ''; }
if(!defined $opt_i) { $opt_i = 4; }
if(!defined $opt_o) { $opt_o = 'dagchainer_rundir'; }
if(!defined $opt_l) { $opt_l = 'cluster_cmds'; }
if(!defined $opt_g) { $opt_g = 'n'; }
if(!defined $opt_p) { $opt_p = 'UGER'; }
if(!defined $opt_q) { $opt_q = 'short'; }
if(!defined $opt_v) { $opt_v = 'n'; }

# Out directory
unless (-d $opt_o) { system("mkdir $opt_o"); }

# Programs
my $DAGCHAINER_PROG = "$Bin/support_scripts/run_DAG_chainer.pl";
my $DAGchainer_to_spans = "$Bin/support_scripts/dagchainer_to_chain_spans.pl";
my $Run_Commands_python = "$Bin/support_scripts/Run_cmds_on_grid.py";
foreach($DAGCHAINER_PROG, $DAGchainer_to_spans, $Run_Commands_python) { die "Cannot find $_ : $!\n" unless (-e $_); }

# Save genomes to compare and those to ignore
my ($data_manager, $genomes_in_repo, $restricted_genomes) = &save_genomes_of_interest($opt_r, $opt_z);

# Save gene_ids from the orthocluster data
my ($orthocluster_to_genes, $genomes_parsed) = &save_gene_ids_from_ortholog_file($opt_c, $restricted_genomes, $genomes_in_repo);

# Process orthocluster results into hit pairs.
my $genome_pair_to_gene_pairs = &process_orthocluster_results_into_hit_pairs($orthocluster_to_genes);

# Split GFF by genome
my $GFF_repo_dictionary = "$opt_r.all.GFF3";
gfffile::split_gff_dictionary_by_species($GFF_repo_dictionary, $opt_r);

# Make dagchainer config file
my $dag_cmds = "-v $opt_v";
my @dagchainer_cmds = &write_dagchainer_conf_file($opt_o, $data_manager, $genomes_parsed, $genome_pair_to_gene_pairs, $dag_cmds);

# Write cmd list
my @cmd_list;
open my $ofh, '>', $opt_l or die "Error, cannot write to $opt_l : $!\n";
foreach my $ind_cmd (@dagchainer_cmds) { 
	print $ofh "$ind_cmd"; 
	push @cmd_list, $ind_cmd;
}
close $ofh;

# Run commands on grid
if($opt_g eq 'y') {
	warn "Running DAGchainer commands on grid...\n";
	my $run_grid_cmd = "$Run_Commands_python --platform $opt_p --queue $opt_q --mem 4 --throttle_nodes 99 --cmds_per_node 1 $opt_l";
	synima::process_cmd($run_grid_cmd);
} 
# Run commands locally
else {
	warn "Running DAGchainer commands...\n";
	foreach(@cmd_list) { system($_); }
}

# Capture the output as a single file
my $dagchainer_output_file = "$opt_r.dagchainer.aligncoords";
my $cmd = "find $opt_o -name \"\*aligncoords\" -exec cat {} + > $dagchainer_output_file";
my $ret = system($cmd);
die "Error, cmd: $cmd died with ret $ret" if ($ret);

# Create spans file
$cmd = "$DAGchainer_to_spans < $dagchainer_output_file > $dagchainer_output_file.spans";
$ret = system($cmd);
die "Error, cmd: $cmd died with ret $ret" if ($ret);

print "Done. Final dagchainer output provided as: $dagchainer_output_file\n";

sub save_genomes_of_interest {
	my ($genomes_interest, $genomes_uninterest) = @_;

	warn "save_genomes_of_interest...\n";
	my %restricted_genomes;
	if(-e $genomes_uninterest) {
		my @genomes = `cat $genomes_uninterest`;
		chomp @genomes;
		%restricted_genomes = map { + $_ => 1 } @genomes;
	}
	my $data_manager = new DataSpecFileParser($genomes_interest);
	my @genomes = $data_manager->get_genome_list();
	my %genomes_in_repo = map { + $_ => 1 } @genomes;
	warn "save_genomes_of_interest: found " . scalar(@genomes) . " genomes\n";
	return ($data_manager, \%genomes_in_repo, \%restricted_genomes);
}

sub save_gene_ids_from_ortholog_file {
	my ($file, $restrict, $include) = @_;

	my (%orthocluster_to_genes, %genomes_parsed);
	open my $fh, '<', $file or die "Error, cannot open file $file : $!\n";
	warn "save_gene_ids_from_ortholog_file: $file...\n";
	CLUSTER: while (my $line=<$fh>) {
		chomp $line;
		next CLUSTER unless ($line =~ m/\w/);
		my @bits = split /\t/, $line;
		my ($cluster_id, $genome_id, $annot_run_name, $trans_id, $gene_id, $locus, $name) = @bits;

		# Ignore these
		next CLUSTER if((%{$restrict}) && (!exists $$restrict{$genome_id}));
		next CLUSTER unless ($$include{$genome_id});

		# Save transcript ids (previously gene_id)
		push (@{$orthocluster_to_genes{$cluster_id}}, { genome => $genome_id, trans_id => "$genome_id:$trans_id"});
		$genomes_parsed{$genome_id}++;
	}
	close $fh;
	warn "save_gene_ids_from_ortholog_file: found " . scalar(keys(%genomes_parsed)) . " genomes\n";
	return (\%orthocluster_to_genes, \%genomes_parsed);
}

sub write_dagchainer_conf_file {
	my ($dagchainer_rundir, $data_manager, $genomes_parsed_href, $genome_pair_to_gene_pairs_href, $dagchainer_commands) = @_;
	my @genomes = sort keys %{$genomes_parsed_href};
	my $number_of_genomes = scalar(@genomes);

	# Build annotation section.
	my @dagchainer_cmds;
	my $dagchainer_cmds_file = "$dagchainer_rundir/dagchainer.cmds";
	open my $cmds_ofh, '>', $dagchainer_cmds_file or die "Cannot open $dagchainer_cmds_file : $!\n";
	warn "Writing DAGChainer config file $dagchainer_cmds_file for $number_of_genomes genomes...\n";	
	for (my $i = 0; $i < $#genomes; $i++) {
		for (my $j = $i + 1; $j <= $#genomes; $j++) {
			my $genome_i = $genomes[$i];
			my $genome_j = $genomes[$j];		
			my $annot_section_text = "";
			my $genome_seq_section_text = "";
			foreach my $genome ($genome_i, $genome_j) {
				my $annotation_file = $data_manager->get_data_dump_filename($genome, "Annotation");
				$annotation_file .= ".synima-parsed.GFF3";
				my $genome_seq_file = $data_manager->get_data_dump_filename($genome, "Genome");
				$annot_section_text .= "$genome = $annotation_file\n";
				$genome_seq_section_text .= "$genome = $genome_seq_file\n";
			}

			my $gene_pairs_aref = $genome_pair_to_gene_pairs_href->{$genome_i}->{$genome_j};
			my $hit_pair_file = "$dagchainer_rundir/${genome_i}_vs_${genome_j}.hit_pairs";
			open my $hp_ofh, '>', $hit_pair_file or die "Cannot open $hit_pair_file : $!\n";
			foreach my $gene_pair_aref (@$gene_pairs_aref) {
				my ($gene_A, $gene_B) = @$gene_pair_aref;
				print $hp_ofh join("\t", $gene_A, $gene_B, "1e-50") . "\n";
			}
			close $hp_ofh;

			# Make DAGchainer config template
			my $dagchainer_template = &get_dagchainer_conf_template($annot_section_text, $genome_seq_section_text, $hit_pair_file, $opt_i);
			my $conf_file = "$dagchainer_rundir/${genome_i}_vs_${genome_j}.dagchainer.conf";
			open my $ofh, '>', $conf_file or die "Error, cannot write to file $conf_file : $!\n";
			print $ofh $dagchainer_template;
			close $ofh;

			# Run DAGchainer
			my $cmd = "$DAGCHAINER_PROG -c $conf_file $dagchainer_commands\n";
			#warn "write_dagchainer_conf_file: $cmd\n";
			print $cmds_ofh $cmd;
			push (@dagchainer_cmds, $cmd);
		}
	}
	close $cmds_ofh;
	return (@dagchainer_cmds);
}

sub process_orthocluster_results_into_hit_pairs {
	my $orthocluster_to_genes = $_[0];
	my %genome_pair_to_gene_pairs;
	foreach my $cluster (keys %{$orthocluster_to_genes}) {
		my @genes = @{$$orthocluster_to_genes{$cluster}};
		for (my $i = 0; $i < $#genes; $i++) {
			for (my $j = $i + 1; $j <= $#genes; $j++) {

				my @pair = ($genes[$i], $genes[$j]);
				@pair = sort {$a->{genome} cmp $b->{genome}} @pair;

				my $genome_A = $pair[0]->{genome};
				my $genome_B = $pair[1]->{genome};
				my $gene_id_A = $pair[0]->{trans_id};
				my $gene_id_B = $pair[1]->{trans_id};
				die "trans_id not saved from clusters. Rerun OrthoMCL_or_RBH_to_summary.pl\n" if((!defined $gene_id_A) || (!defined $gene_id_B));
				push (@{$genome_pair_to_gene_pairs{$genome_A}->{$genome_B}}, [$gene_id_A, $gene_id_B]);
			}
		}
	}
	warn "process_orthocluster_results_into_hit_pairs: found " . scalar(keys(%genome_pair_to_gene_pairs)) . " genomes\n";
	return \%genome_pair_to_gene_pairs;
}

sub get_dagchainer_conf_template {
	my ($annotation_file_def, $genome_seq_def, $hit_pairs_def, $min_pairs) = @_;	
	my $conf_template = "#-----------------------------------------------------------------------
	;; configuration file for dagchainer

[GeneAnnotations]

$annotation_file_def

[GenomeSequences]

$genome_seq_def

[MatchPairs]

# only first three fields of a tab-delimited file are examined, expecting:
#  accA accB E-value
# compatible with NCBI-blast (blastall)  -m 8 output format.

Data = $hit_pairs_def

[NoiseFilter]
BEST_MATCH_AGGREGATE_DIST = 5

[Orthologs]
# Data = all_orthomcl.out

# Section is optional

[Parameters]

# MODE can correspond to relative gene position or actual genome coordinate
# with respective values:  RELATIVE_POSITION  |   GENOME_COORDINATE
MODE = RELATIVE_POSITION

# gap open and extend penalties
GAP_OPEN = 0
GAP_EXTEND = -1

## size of a single gap
GAP_LENGTH = 1

MAX_MATCH_SCORE = 50

# comment out the line below to enforce a constant match score 
# instead of the min(-log(Evalue), MAX_MATCH_SCORE) value.
CONSTANT_MATCH_SCORE = 3

# maximum E-value 
MAX_EVALUE = 10

# maximum distance allowed between two neighboring syntenic genes in a single block
MAX_DIST_BETWEEN_SYN_PAIRS = 5


# minimum alignment score of the highest scoring block
# by default, this is set dynamically to: MIN_ALIGN_LEN * 2.5 * -GAP_PENALTY
# MIN_ALIGNMENT_SCORE = 

# minimum number of aligned gene pairs within a single block
MIN_ALIGNED_PAIRS = $min_pairs

# Include self-molecule comparisons.  Turn on if looking for segmental genome duplications
INCLUDE_SELF_COMPARISONS = FALSE

VERBOSE = FALSE

; #------------------------------------------------------------------------------------------";
	return($conf_template);
}
