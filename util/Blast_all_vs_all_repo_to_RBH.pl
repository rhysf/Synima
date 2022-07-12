#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use read_FASTA;
use Synima;
use DataSpecFileParser;
use RBH_manager;
use File::Basename;
use Cwd;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <repo_spec>\n
Optional: -t Type (PEP/CDS) [PEP]
	  -o Out directory [RBH_outdir]
	  -g Run commands on the grid (y/n) [n]
	  -p Platform (UGER, LSF, GridEngine) [UGER]
	  -q Queue name [short]\n";
our($opt_r, $opt_t, $opt_o, $opt_g, $opt_p, $opt_q);
getopt('rtogpq');
die $usage unless($opt_r);
if(!defined $opt_t) { $opt_t = 'PEP'; }
if(!defined $opt_o) { $opt_o = 'RBH_outdir'; }
if(!defined $opt_g) { $opt_g = 'n'; }
if(!defined $opt_p) { $opt_p = 'UGER'; }
if(!defined $opt_q) { $opt_q = 'short'; }
die "-g is not n, N, y or Y: $opt_g\n" if($opt_g !~ m/^(n|N|y|Y)$/);
die "-t is not PEP or CDS: $opt_t\n" if($opt_t !~ m/^(PEP|CDS)$/);
die "-p is not UGER, LSF or GridEngine: $opt_p\n" if($opt_p !~ m/^(UGER|LSF|GridEngine)$/);

# Dependencies
my $uname = `uname`;
chomp $uname;
my $slclust = "$Bin/support_scripts/slclust.$uname";
foreach($slclust) { die "Cannot find $_ : $!\n" unless(-e $_); }

# Write RBH files
my $RBH_commands2 = &write_RBH_files($opt_r, $opt_t, $opt_o);
die "No commands made in $RBH_commands2. Check input files\n" if((-s $RBH_commands2) eq 0);

# Process
my $cmd = "cut -f1,2 $RBH_commands2 > $opt_o/$opt_t.RBH.pairs";
synima::process_cmd($cmd);

# Cluster using PASA C++ program
$cmd = "$slclust < $opt_o/$opt_t.RBH.pairs > $opt_o/$opt_t.RBH.slclust";
synima::process_cmd($cmd);

# Add in paralogs
&add_in_paralogs_to_orthoclusters($opt_r, "$opt_o/$opt_t.RBH.slclust", $opt_t, "$opt_o/$opt_t.RBH.OrthoClusters");
warn "Finished.\n";

sub write_RBH_files {
	my ($repo_spec, $type, $output) = @_; 
	warn "write_RBH_files: Printing to $output\n";

	# Make output directory
	`mkdir $output` if(! -d $output);
	my $RBH_commands_cat = "$output/$type.RBH";
	if(-e $RBH_commands_cat) { 
		warn "Warning, $RBH_commands_cat already exists. Not overwriting.\n";
		return $RBH_commands_cat;
	}

	# Genomes in lexical order from repo spec
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();
	@genomes = sort @genomes; 
	die "No genomes found in $repo_spec\n" if($#genomes eq 0);

	# make sure each protein file is blastable
	for (my $i = 0; $i <= $#genomes; $i++) {
		my $genomeA = $genomes[$i];
		GENOMES2: for (my $j = $i; $j <= $#genomes; $j++) {
			my $genomeB = $genomes[$j];
			my $genomeA_vs_genomeB_blast_file = &RBH_manager::get_blast_output_file($data_manager, $genomeA, $genomeB, $type);
			my $genomeB_vs_genomeA_blast_file = &RBH_manager::get_blast_output_file($data_manager, $genomeB, $genomeA, $type);
			die "Error, no required blast output file: $genomeA_vs_genomeB_blast_file" if (! -s $genomeA_vs_genomeB_blast_file);
			die "Error, no required blast output file: $genomeB_vs_genomeA_blast_file" if (! -s $genomeB_vs_genomeA_blast_file);

			# Print Reciprical best matches
			my $genomes_RBH_file = &RBH_manager::get_RBH_output_file($data_manager, $genomeA, $genomeB, $type);
			if (-s $genomes_RBH_file) { print STDERR "Warning, $genomes_RBH_file already exists, not replacing it.\n"; }
			else { &m8_recip_best_match($genomeA_vs_genomeB_blast_file, $genomeB_vs_genomeA_blast_file, $genomes_RBH_file); }

			# Concatenate
			open my $fh, '<', $genomes_RBH_file or die "Cannot open file $genomes_RBH_file: $!\n";
			open my $ofh, '>>', $RBH_commands_cat or die "Cannot open $RBH_commands_cat: $!\n";
			while (my $line=<$fh>) { print $ofh $line; }
		   	close $fh;
			close $ofh;
		}
	}
	return $RBH_commands_cat;
}

sub m8_recip_best_match {
	my ($file1, $file2, $out) = @_;

	# Save the best matches in structure
	my %best_matches;
	foreach my $m8_file ($file1, $file2) {
		open my $fh, '<', $m8_file or die "Cannot open $m8_file : $!\n";
		LINES: while (my $line=<$fh>) {
			chomp $line;
			my @x = split /\t/, $line;
	   		die "$m8_file misformed ($line). Run from start\n" unless(defined $x[11]);
			my ($acc_A, $acc_B, $bit) = ($x[0], $x[1], $x[11]);
			next LINES if($acc_A eq $acc_B);

			my $best_score_A = $best_matches{$acc_A};
			if ((! $best_score_A) || ($bit > $best_score_A->{bit})) {
		   		$best_matches{$acc_A} = { 
					bit => $bit,
					hits => [ {line => $line, acc => $acc_B}, ],
				   };
			}

			# Tied score
			elsif ($bit == $best_score_A->{bit}) {
		   		push (@{$best_matches{$acc_A}->{hits}}, { 
						line => $line,
						acc => $acc_B,
					});
			}
	   	}
	}

	# Print
	open my $fh, '>', $out or die "Cannot open $out: $!\n";
	my @accs = keys %best_matches;
	foreach my $acc (@accs) {
		if(my $struct = $best_matches{$acc}) {
			my @hits = @{$struct->{hits}};
			foreach my $hit (@hits) {
				my $other_acc = $hit->{acc};
				if (my $recip_hit = &find_recip_hit($other_acc, $acc, \%best_matches)) {
					print $fh $hit->{line} . "\t<AND>\t" . $recip_hit->{line} . "\n";
				}
			}
			delete $best_matches{$acc};
		}
	}
	close $fh;
	return 1;
}

sub find_recip_hit {
	my ($other_acc, $acc, $best_matches_href) = @_;

	my $struct = $best_matches_href->{$other_acc};
	unless (ref $struct) {
		return undef;
	}

	my @hits = @{$struct->{hits}};
	foreach my $hit (@hits) {
		if ($hit->{acc} eq $acc) {
			return($hit);
		}
	}
	return undef;
}

sub add_in_paralogs_to_orthoclusters {
	my ($repo_spec, $slclust_file, $type, $output) = @_; 
	my $data_manager = new DataSpecFileParser($repo_spec);

	warn "--getting gene annotations\n";
	my $annot_id_to_annotations = &get_gene_annotations($data_manager);
	
	warn "--parsing clusters from $slclust_file\n";
	my $cluster_id_to_clusters = &parse_clusters($slclust_file);

	warn "--mapping genes to their cluster IDs.\n";
	my $gene_to_cluster = &map_gene_to_cluster_id($cluster_id_to_clusters);

	warn "--getting top blast score per orthologous gene\n";
	my $gene_to_top_ortho_blast_score = &get_top_ortho_blast_score($data_manager, $type);

	warn "--capturing InParalogs\n";
	my $cluster_id_to_InParalogs = &get_InParalogs($data_manager, $gene_to_top_ortho_blast_score, $gene_to_cluster, $type);

	## report final results including InParalogs
	open my $ofh, '>', $output or die "Cannot open output file $output : $!\n";
	foreach my $ortho_cluster_id (sort {$a<=>$b} keys %{$cluster_id_to_clusters}) {

		# Orthologs
		my @orthologs = @{$$cluster_id_to_clusters{$ortho_cluster_id}};
		foreach my $ortholog (@orthologs) {
			my $gene_struct = $$annot_id_to_annotations{$ortholog} or die "Error, no annotation struct for [$ortholog] : $!\n";
			my $genome = $gene_struct->{genome};
			my $gene_id = $gene_struct->{gene_id};
			#my $transcript_id = $gene_struct->{trans_id};
			#my $locus_id = $gene_struct->{locus};
			my $name = $gene_struct->{name};
			print $ofh "$ortho_cluster_id\tOrtho\t$genome\t$gene_id\t$gene_id\t$gene_id\t$name\n";
		}

		# Paralogs
		if (exists $$cluster_id_to_InParalogs{$ortho_cluster_id}) {
			my @in_paralogs = @{$$cluster_id_to_InParalogs{$ortho_cluster_id}};
			foreach my $in_paralog (@in_paralogs) {
				my $gene_struct = $$annot_id_to_annotations{$in_paralog} or die "Error, no annotation struct for InParalog: [$in_paralog] :$!\n";
				my $genome = $gene_struct->{genome};
				my $gene_id = $gene_struct->{gene_id};
				#my $transcript_id = $gene_struct->{trans_id};
				#my $locus_id = $gene_struct->{locus};
				my $name = $gene_struct->{name};
				print $ofh "$ortho_cluster_id\tInPara\t$genome\t$gene_id\t$gene_id\t$gene_id\t$name\n";
			}
		}

		# spacer between groups
		print $ofh "\n"; 
	}
	close $ofh;
	return 1;
}

sub get_gene_annotations {
	my ($data_manager) = @_;

	my @genomes = $data_manager->get_genome_list();
	my %annot_id_to_annotations;
	foreach my $genome (@genomes) {

		# Save sequences
		my $peptide_file = $data_manager->get_data_dump_filename($genome, 'PEP');
		my $sequences = fastafile::fasta_to_struct($peptide_file);
		foreach my $id(keys %{$$sequences{'seq'}}) {
			my $struct = { 
				gene_id => $id,
				trans_id => $id,
				genome => $genome,
				locus => $id,
				name => $$sequences{'desc'}{$id},
			};
			$annot_id_to_annotations{$id} = $struct;
		}
	}
	return \%annot_id_to_annotations;
}

sub parse_clusters {
	my ($slclust_file) = @_;

	my $counter = 0;
	my %cluster_id_to_clusters;
	open (my $fh, $slclust_file) or die "Error, cannot open file $slclust_file";
	while(my $line=<$fh>) {
		chomp $line;
		$counter++;
		my @x = split /\s+/, $line;
		push (@{$cluster_id_to_clusters{$counter}}, @x);
	}
	close $fh;
	return \%cluster_id_to_clusters;
}

sub map_gene_to_cluster_id {
	my ($cluster_id_to_clusters_href) = @_;
	my %gene_to_cluster_id;
	foreach my $cluster_id (keys %$cluster_id_to_clusters_href) {
		my @genes = @{$cluster_id_to_clusters_href->{$cluster_id}};
		foreach my $gene (@genes) { $gene_to_cluster_id{$gene} = $cluster_id; }
	}
	return \%gene_to_cluster_id;
}

sub get_top_ortho_blast_score {
	my ($data_manager, $type) = @_;
	
	my @genomes = $data_manager->get_genome_list();
	my %id_to_top_blast_score;
	for (my $i = 0; $i <$#genomes; $i++) {
		my $genomeA = $genomes[$i];
		for (my $j = ($i + 1); $j <= $#genomes; $j++) {
			my $genomeB = $genomes[$j];
			my $genomes_RBH_file = &RBH_manager::get_RBH_output_file($data_manager, $genomeA, $genomeB, $type);

			die "Error, $genomes_RBH_file does not exist! : $!\n" if (! -s $genomes_RBH_file);
			open my $fh, '<', $genomes_RBH_file or die "Error, cannot open file $genomes_RBH_file : $!\n";
			while (my $line=<$fh>) {
				chomp $line;
				my @x = split /\t/, $line;
				my ($accA, $accB, $bit_score) = ($x[0], $x[1], $x[11]);

				# Save BLAST bit scores, or overwrite smaller scores
				if ((! exists $id_to_top_blast_score{$accA}) || ($id_to_top_blast_score{$accA} < $bit_score)) {
					$id_to_top_blast_score{$accA} = $bit_score;
				}
				if ((! exists $id_to_top_blast_score{$accB}) || ($id_to_top_blast_score{$accB} < $bit_score)) {
		        		$id_to_top_blast_score{$accB} = $bit_score;
				}
			}
			close $fh;
		}
	}
	return \%id_to_top_blast_score;
}

sub get_InParalogs {
	my ($data_manager, $gene_to_top_ortho_blast_score_href, $gene_to_cluster_href, $type) = @_;

	my %inParalog_to_ortho_group;
	my @genomes = $data_manager->get_genome_list();
	foreach my $genome (@genomes) {
		my $self_blast_file = &RBH_manager::get_blast_output_file($data_manager, $genome, $genome, $type);
		open my $fh, '<', $self_blast_file or die "Error, cannot open file $self_blast_file : $!\n";
		while (my $line=<$fh>) {
			chomp $line;
			my @x = split/\t/, $line;
			my ($accA, $accB, $bit_score) = ($x[0], $x[1], $x[11]);

			# Ignore self hits
			next if ($accA eq $accB);

		 	# want a case where one is in an ortho group, the other is not, and the match between them is better than the best ortho match.
		 	# want one or the other, not both or neither
		 	next unless ( (exists ($gene_to_top_ortho_blast_score_href->{$accA})) xor (exists ($gene_to_top_ortho_blast_score_href->{$accB})));

		 	# Swap them so that accA is the one in the ortho cluster
			if (exists ($gene_to_top_ortho_blast_score_href->{$accB})) { ($accA, $accB) = ($accB, $accA); }
		 
			# Got one. B matches A better than A matches its best-matching ortholog
			if ($bit_score > $gene_to_top_ortho_blast_score_href->{$accA}) {
				my $ortho_cluster = $gene_to_cluster_href->{$accA} or die "Error, no cluster found for acc: $accA\n";
				if (my $struct = $inParalog_to_ortho_group{$accB}) {
					if (($struct->{cluster} != $ortho_cluster) && ($struct->{score} < $bit_score)) {
						print STDERR "Warning, $accB already mapped to orthogroup: $struct->{cluster} with bit score: $struct->{score}, but scores higher when mapped to cluster $ortho_cluster with score $bit_score. Moving it.\n";
						$inParalog_to_ortho_group{$accB} = {cluster => $ortho_cluster, score => $bit_score, };
					}
				} else {	
					$inParalog_to_ortho_group{$accB} = { cluster => $ortho_cluster, score => $bit_score, };
				}
			}
		}
		close $fh;
	}

	# Transform the data structure to this: ortho_id => [geneA, geneB, ...]
	my %cluster_id_to_para_list;
	foreach my $acc (keys %inParalog_to_ortho_group) {
		my $ortho_group_struct = $inParalog_to_ortho_group{$acc};
		my $ortho_group = $ortho_group_struct->{cluster};
		push (@{$cluster_id_to_para_list{$ortho_group}}, $acc);
	}
	return \%cluster_id_to_para_list;
}
