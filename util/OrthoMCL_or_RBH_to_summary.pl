#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use read_FASTA;
use Ort;
use Carp qw (cluck croak confess);
use Cwd;
use Synima;

### rfarrer@broadinstitute.org

# Opening statements
my $usage = "Usage: perl $0 -o <ortholog file (E.g. PEP.RBH.OrthoClusters, all_orthomcl.out)>\n
Optional -t Type of clustering (OMCL, RBH) [OMCL]
         -d Outdir from Blast_all_vs_all_repo_to_OrthoMCL/RBH.pl [OMCL_outdir]
	 -r Repo Spec [./Repo_spec.txt]
	 -p Repo Spec Peptide file [./Repo_spec.txt.all.pep]\n";
our($opt_t, $opt_d, $opt_r, $opt_p, $opt_o);
getopt('tdrpo');
if(!defined $opt_t) { $opt_t = 'OMCL'; }
if(!defined $opt_d) { $opt_d = 'OMCL_outdir'; }
if(!defined $opt_r) { $opt_r = './Repo_spec.txt'; }
if(!defined $opt_p) { $opt_p = './Repo_spec.txt.all.pep'; }
die $usage unless($opt_o);
die "$opt_o not recognised as a file\n" if(! -f $opt_o);
die "-t not equal to OMCL or RBH: $opt_t\n" unless(($opt_t eq 'OMCL') || ($opt_t eq 'RBH'));
my $genome_codes;
if($opt_t eq 'OMCL') { 
	die "Cannot find outdir $opt_d : $!\n" unless(-d $opt_d);
	$genome_codes = "$opt_d/for_omcl.genome_codes";
	die "Cannot find $genome_codes : $!\n" unless(-e $genome_codes);
} else { $genome_codes = 'NA'; }

# Parse repo spec into struct
warn "Parsing repo spec...\n";
my $data_manager = new DataSpecFileParser($opt_r);

# Write all outputs to the orthocluster out directory in the repo
my $ortholog_analysis_name = 'GENE_CLUSTERS_SUMMARIES';
my $ortholog_analysis_dir = ($ortholog_analysis_name . '.' . $opt_t);
my $repo_dir_name = $data_manager->get_repo_dir();
my $orthocluster_out_dir = "$repo_dir_name/$ortholog_analysis_dir";
if (! -d $orthocluster_out_dir) { mkdir ($orthocluster_out_dir) or die "Error, cannot mkdir $orthocluster_out_dir : $!\n"; }
my $functional_annotation_map = ($orthocluster_out_dir . '/transcript_annotation.map');

# Dump the ortholog report.
warn "Reading ORTs file...\n";
my $orts = Ort::new($opt_o, 'transcript', $opt_t, $genome_codes);
$orts->read();
fastafile::create_transcript_functinal_annotation_map($opt_p, $functional_annotation_map);
$orts->add_annotation($functional_annotation_map);

# Write new ortholog file
my $clusters_outfile = "$orthocluster_out_dir/$ortholog_analysis_name.clusters";
warn "Writing new ORTs file $clusters_outfile...\n";	
open my $ofh, '>', $clusters_outfile or die "Unable to write on $clusters_outfile file : $!\n";
foreach my $cluster_num ( $orts->get_cluster_ids() ) {
	print $ofh $orts->get_cluster_desc_repo_format($cluster_num);
	print $ofh "\n";
}
close $ofh;

# Save protein struct
my $all_proteins = fastafile::fasta_broad_format_to_struct($opt_p); 

# Dump unique gene report:
warn "Generating uniques file...\n";
my $num_uniq_transcripts = 1;
my $unique_genes_file = "$orthocluster_out_dir/$ortholog_analysis_name.unique";
open $ofh, '>', $unique_genes_file or die "Unable to write on $unique_genes_file file : $!\n";
foreach my $curr_transcript ( keys %{$all_proteins}) {
	if( not $orts->has_at_least_one_ort( $curr_transcript ) ) {
		print $ofh "uniq_" . $num_uniq_transcripts . "\t" .
			$$all_proteins{$curr_transcript}{genome} . "\t" .
			"Ortho" . "\t" .
			$$all_proteins{$curr_transcript}{transcript_id} . "\t" .
			$$all_proteins{$curr_transcript}{gene_id} . "\t" .
			$$all_proteins{$curr_transcript}{locus_name} . "\t" .
			$$all_proteins{$curr_transcript}{func_annot} . "\n";
		$num_uniq_transcripts++;				 
	}
}
close $ofh;

# Combine 'cos it's useful
`cat $orthocluster_out_dir/$ortholog_analysis_name.clusters $orthocluster_out_dir/$ortholog_analysis_name.unique >> $orthocluster_out_dir/$ortholog_analysis_name.clusters_and_uniques`;

# Summarize counts of orthologous genes per gene cluster per genome
warn "Summarizing counts of orthologous genes per gene cluster per genome...\n";
my $cluster_counts_file = "$orthocluster_out_dir/$ortholog_analysis_name.cluster_dist_per_genome.txt";
&count_clustered_genes_per_genome($clusters_outfile, $unique_genes_file, 'transcript', $cluster_counts_file);

warn "Plotting cluster results...\n";	
&omcl_matrix_to_count_stats($cluster_counts_file, "$cluster_counts_file.summary", $orthocluster_out_dir);

warn "Plotting cluster heatmap...\n";
&omcl_matrix_to_heatmap($cluster_counts_file, $orthocluster_out_dir);

warn "Generating reference genome views of cluster content...\n";
&omcl_matrix_by_ref_genome_view($opt_r, $cluster_counts_file, $orthocluster_out_dir);

# usage: $0 gene_clusters.txt unique_genes.txt < locus identifer = 'transcript' or 'gene' >
sub count_clustered_genes_per_genome {
	my ($clusters, $unique_genes, $locus_identifier, $outfile) = @_;

	# Go through clusters and unique genes saving info
	my (%cluster_to_genome_count, %genomes, %cluster_to_name_count, %seen, %genome_to_gene_count);
	FILES: foreach my $file ($clusters, $unique_genes) {
		open my $fh, '<', $file or die "Error, cannot open file $file : $!\n";
		LINES: while (my $line=<$fh>) {
			chomp $line;
			next LINES unless ($line =~ m/\w/);
			my @x = split /\t/, $line;
			my ($cluster_id, $genome, $name) = ($x[0], $x[1], $x[6]);
			die "Bad format. Name not defined in $line in $file. Re-make\n" if(!defined $name);
			my $id;
			if($locus_identifier eq 'gene') { $id = $x[4]; }
			else { $id = $x[3]; }

			#die "Error, $locus_identifier $id shows up multiple times in $file\n" if($seen{$id});
			#$seen{$id}++;

			$genomes{$genome}++;
			$cluster_to_genome_count{$cluster_id}{$genome}++;
			$cluster_to_name_count{$cluster_id}{$name}++;
			$genome_to_gene_count{$genome}++;
		}
		close $fh;
	}
	my @organisms = sort keys %genomes;

	# Outfile
	open my $ofh, '>', $outfile or die "Cannot open $outfile: $!\n";

	# count header
	my $count_header = "#";
	foreach my $genome (keys %genome_to_gene_count) {
		my $count = $genome_to_gene_count{$genome};
		$count_header .= "$genome=$count\t";
	}
	chomp $count_header;
	print $ofh "$count_header\n";

	# header
	print $ofh "#cluster_id\tname\t" . join("\t", @organisms) . "\n";

	foreach my $cluster (sort keys %cluster_to_genome_count) {
		print $ofh "$cluster";
		my $names_href = $cluster_to_name_count{$cluster};
	
		my @names = sort {$names_href->{$a}<=>$names_href->{$b}} keys %$names_href;
		my $most_popular_name = pop @names;
		$most_popular_name =~ s/\W/_/g; # non-word chars cause problems with some matrix analysis tools
		print $ofh "\t$most_popular_name";

		foreach my $organism (@organisms) {
			my $count = $cluster_to_genome_count{$cluster}->{$organism} || 0;
			print $ofh "\t$count";
		}
		print $ofh "\n";
	}
	close $ofh;
	return 1;
}

sub omcl_matrix_to_count_stats {
	my ($cluster_counts_file, $output_file, $output_dir) = @_;
	my %genomes;
	my $genomes_href = \%genomes;
	my $genomes_list_file; # genomes_list_file.txt

	if ($genomes_list_file) {
		my @genome_list = `cat $genomes_list_file`;
		foreach my $genome (@genome_list) {
			$genome =~ s/\s//g;
			$genomes{$genome} = 1 if $genome =~ /\w/;
		}
	}

	my %field_to_genome;
	open my $fh, '<', $cluster_counts_file or die "Error, cannot open file $cluster_counts_file : $!\n";
	{
		my $trash = <$fh>; ####
		my $header = <$fh>;
		chomp $header;
		
		my @fields = split(/\t/, $header);
		# first two fields are gene_id and name
		shift @fields;
		shift @fields;
		
		for (my $i = 0; $i <= $#fields; $i++) {
		 	my $genome = $fields[$i];
		 	if ( (! %$genomes_href) || exists $genomes_href->{$genome}) {
		 		$field_to_genome{$i} = $fields[$i];
		 	} else { warn "-skipping genome: $genome\n"; }
		}
	}
	
	my %genome_to_class_count;
	my %genome_type_to_cluster_ids;

	while (<$fh>) {
		chomp;
		my ($cluster_id, $name, @genome_counts) = split(/\t/);
		
		my @counts;
		foreach my $field (keys %field_to_genome) {
		my $count = $genome_counts[$field];
		my $genome = $field_to_genome{$field};
		
		push (@counts, { genome => $genome,
				count => $count,
		 });
		}

		my @counts_gt_zero = grep { $_->{count} > 0 } @counts;

		my $classification = "aux";
		if (scalar @counts_gt_zero == scalar(@counts)) { $classification = "core"; }
		elsif (scalar @counts_gt_zero == 1) { $classification = "unique"; }
		else { }
		
		# print "$classification: @counts_gt_zero\n";
		foreach my $count_struct (@counts) {
			my $count = $count_struct->{count};
			my $genome = $count_struct->{genome};
			if ($count) {
				$genome_to_class_count{$genome}->{$classification}+=$count;
				$genome_type_to_cluster_ids{$genome}->{$classification}->{$cluster_id}++;
			}
		}
		
	}
	close $fh;

	my @genomes = sort keys %genome_to_class_count;
	
	my @core_counts;
	my @aux_counts;
	my @uniq_counts;

	# Summary output
	open OUT, ">$output_file" or die "Cannot open $output_file\n";
	print OUT join("\t", "#genome", "core", "aux", "unique") . "\n";
	foreach my $genome (@genomes) {
		my $core_count = $genome_to_class_count{$genome}->{'core'} || 0;
		my $aux_count = $genome_to_class_count{$genome}->{'aux'} || 0;
		my $uniq_count = $genome_to_class_count{$genome}->{'unique'} || 0;
				
		push (@core_counts, $core_count);
		push (@aux_counts, $aux_count);
		push (@uniq_counts, $uniq_count);
	
		print OUT join("\t", $genome, $core_count, $aux_count, $uniq_count) . "\n";
	}

	print OUT "// Classifications:\n";
	foreach my $genome (keys %genome_type_to_cluster_ids) {
		my $type_href = $genome_type_to_cluster_ids{$genome};
		foreach my $type (sort keys %$type_href) {
			my $gene_ids_href = $type_href->{$type};
			foreach my $gene_id (keys %$gene_ids_href) {
				print OUT "$genome\t$type\t$gene_id\n";
			}
			}
		}
		close OUT;

	# Rscript should go in the ortholog folder. Not CWD
	my $rscript = ($output_dir . "/plot_cluster_dist.Rscript");
	open my $ofh, '>', $rscript or die "Error, cannot write to file $rscript : $!\n";
	print $ofh "core_counts = c(" . join(",", @core_counts) . ")\n";
	print $ofh "aux_counts = c(" . join(",", @aux_counts) . ")\n";
	print $ofh "uniq_counts = c(" . join(",", @uniq_counts) . ")\n";
	print $ofh "data = cbind(core_counts, aux_counts, uniq_counts)\n";
	print $ofh "data = t(data)\n";
	print $ofh "genomeNames = c(\'" . join("','", @genomes) . "\')\n";
	print $ofh "postscript(file=\"$output_dir/cluster_dist.eps\")\n";
	print $ofh "par(mar = c(10, 4, 4, 2))\n";
	print $ofh "barplot(data, main=\"Dist. for Ortholog Clusters\", names=genomeNames, col=c('green', 'blue', 'yellow'), las=2, cex.names=0.8)\n";
	print $ofh "dev.off()\n";
	close $ofh;
	
	my $cmd = "R --vanilla -q 1>&2  < $rscript";
	synima::process_cmd($cmd);
	return;
}

sub omcl_matrix_to_heatmap {
	my ($ortho_file, $outdir) = @_;
	my $rscript = "$outdir/write_ortho_heatmap.Rscript";
	open my $ofh, '>', $rscript or die "Error, cannot write to $rscript : $!\n";
	print $ofh "data = read.table(\"$ortho_file\", com=\'\', sep=\'\\t\', skip=1, header=T)\n";
	print $ofh "row.names(data) = paste(data[,1], data[,2])\n";
	print $ofh "data = data[,3:length(data[1,])]\n";
	print $ofh "data_matrix = data.matrix(data)\n";
	print $ofh "postscript(\"$ortho_file.heatmap.eps\")\n";
	print $ofh "par\(yaxt=\'n\'\) \# turn off gene labels\n";
	print $ofh "data_heatmap = heatmap(data_matrix, Rowv=NA, scale=\"none\", col=c('black', 'green', 'red'), zlim=c(0,2), margins=c(12,1))\n";
	print $ofh "dev.off()\n";
	close $ofh;
	my $cmd = "R --vanilla -q 1>&2  < $rscript";
	synima::process_cmd($cmd);
}

sub omcl_matrix_by_ref_genome_view {
	my ($repo_file, $clusters_file, $outdir_main) = @_;
	my $data_manager = new DataSpecFileParser($repo_file);

	# Make out directory
	my $outdir = "$outdir_main/clusters_arranged_by_ref_genome";
	if (! -d $outdir) { mkdir $outdir or die "Error, cannot mkdir $outdir : $!\n"; }

	# gene_to_omcl =  gene => cluster_id
	# omcl_to_genome_gene = cluster_id => genome => { geneA => 1, geneB=> 1 }
	my ($gene_to_omcl, $omcl_to_genome_gene) = &parse_omcl_clusters($clusters_file);

	# Save GFF
	my @genomes = $data_manager->get_genome_list();
	my %genome_to_annots;
	foreach my $genome (@genomes) {
		my $annot_file = $data_manager->get_data_dump_filename($genome, 'Annotation');
		warn "-processing $annot_file\n";
		&parse_annot($genome, $annot_file, \%genome_to_annots);
	}

	my @ref_genomes = keys %genome_to_annots;
	foreach my $refGenomeName (@ref_genomes) {

		my $outfile = "$outdir/$refGenomeName.cluster_content_compare";
		open my $ofh, '>', $outfile or die "Error, cannot write to file $outfile : $!\n";
		select $ofh;
		
		warn "// writing to $outfile\n";
		
		## generate reference-based matrix:
		my @other_genomes = grep { $_ ne $refGenomeName } keys %genome_to_annots;
		my @ref_genes = sort {$a->{scaff} cmp $b->{scaff} || $a->{lend} <=> $b->{lend} } @{$genome_to_annots{$refGenomeName}};

		print join("\t", "#annot", $refGenomeName, @other_genomes) . "\n";
		foreach my $ref_gene (@ref_genes) {
			
			my $annot = $ref_gene->{scaff};
			if ($ref_gene->{locus}) { $annot .= "-" . $ref_gene->{locus}; }
			$annot .= "-" . $ref_gene->{name};
			$annot =~ s/\s+/_/g;
			my $ref_gene_id = $ref_gene->{gene_id};
			my $omcl_id = $$gene_to_omcl{$ref_gene_id};
			
			if ($omcl_id) {
				my $other_genomes_clusters_href = $$omcl_to_genome_gene{$omcl_id};
				
				my $ref_genome_count = 0;
				my $ref_genes_href = $other_genomes_clusters_href->{$refGenomeName};
				if ($ref_genes_href) {
					$ref_genome_count = scalar (keys %$ref_genes_href);
				}
				
				print "$annot-omcl($omcl_id)\t" . $ref_genome_count;
				
				foreach my $other_genome (@other_genomes) {
					my $count = 0;
					if (my $other_genes_href = $other_genomes_clusters_href->{$other_genome}) {
						$count = scalar(keys %$other_genes_href);
					}
					print "\t$count";
				}
				print "\n";
			}
			else {
				#print "no omcl_id for gene: $ref_gene_id\n";
				print "$annot\t1";
				# not in an omcl
				foreach my $other_genome (@other_genomes) { print "\t0"; }
				print "\n";
			}
		}
		close $ofh;
	}
	return 1;
}

sub parse_omcl_clusters {
	my ($clusters_file) = $_[0];
	my (%gene_to_omcl, %omcl_to_genome_gene);
	open my $fh, '<', $clusters_file or die "Error, cannot open file $clusters_file : $!\n";
	while (my $line=<$fh>) {
		next if($line =~ /^#/);
		next unless ($line =~/\w/);
		my @x = split /\t/, $line;
		my ($cluster_id, $genome_name, $analysisRunName, $trans_id, $gene_id, $locus, $com_name) = @x;
		$gene_to_omcl{$gene_id} = $cluster_id;
		$omcl_to_genome_gene{$cluster_id}{$genome_name}{$gene_id} = 1;
	}
	close $fh;
	return (\%gene_to_omcl, \%omcl_to_genome_gene);
}

sub parse_annot {
	my ($genome, $annot_file, $genome_to_annots_href) = @_;

	warn "Saving info from $annot_file...\n";
	open (my $fh, $annot_file);
	while (my $line=<$fh>) {
		next unless ($line =~ m/\w/); 
		my @x = split /\t/, $line;
		my ($scaff, $source, $type, $lend, $rend, $x1, $orient, $x2, $info) = @x;
		next unless($type eq "gene");

		# Save gene id, Name and alias from GFF3
		my ($gene_id)  = ( $info =~ /ID=([\w\W]+?)[;\n]/ );
		die "Error, cannot parse ID info from $info:\n$line\n" if(!defined $gene_id);
		my ($name) = ( $info =~ /Name=([\w\W]+?)[;\n]/ );
		die "Error, cannot parse name info from $info $name:\n$line\n" if(!defined $name);
		my ($alias) = ( $info =~ /Alias=([\w\W]+?)[;\n]/ );
		my $locus;
		if (defined $alias) {	
			$locus = $gene_id;
			$gene_id = $alias; # replace with calhoun gene identifier
		}
		my $gene_info = { scaff => $scaff,
			   	 lend => $lend,
			   	 rend => $rend,
			   	 gene_id => $gene_id,
			   	 name => $name,
			   	 orient => $orient,
			   	 locus => $locus || "",
			};
		push (@{$genome_to_annots_href->{$genome}}, $gene_info);
	}
	close $fh;
	return;
}
