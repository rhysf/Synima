package gfffile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);

### rfarrer@broadinstitute.org

sub gff_to_contig_parent_to_cds_hash {
	my ($input) = @_;
	my (%hash_info, %hash_strand);
	warn "gff_to_contig_parent_to_cds_hash: saving all from $input...\n";
	open my $fh, '<', $input or die "Cannot open $input: $!\n";
	GFF: while(my $line = <$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		next GFF if ((@bits < 8) || ($line =~ m/^\#/));
		my ($contig, $source, $feature, $start, $stop, $score, $strand, $frame, $gene) = @bits;
		my $CDS = "$start-$stop";

		# Save (avoid leaving a space at the end)
		if(!defined $hash_info{$contig}{$gene}) { $hash_info{$contig}{$gene} = $CDS; }
		else { $hash_info{$contig}{$gene} .= " $CDS"; }
		$hash_strand{$gene} = $strand;
	}
	close $fh;
	return (\%hash_info, \%hash_strand);
}

sub combine_all_gff3_files_in_repo {
	my ($repo_spec, $output, $feature, $desc_seperator, $desc_column, $desc_replace) = @_;
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();

	# Unique transcript ids needed
	my %unique_transcript_ids;
	my %new_transcript_ids;

	# Combine all GFFs
	warn "combine_all_gff3_files_in_repo: printing to $output\n";
	die "$output already exists. Delete then re-run\n" if(-e $output);
	open my $ofh, '>>', $output or die "Error, cannot open file $output\n";
	foreach my $genome (@genomes) {

		# open GFF for each genome:
		my %ids_in_one_genome;
		my $count = 0;
		my $change_id_count = 0;
		my $duplicate_id_count = 0;
		my $annot_gff3 = $data_manager->get_data_dump_filename($genome, "Annotation");
		warn "combine_all_gff3_files_in_repo: opening $annot_gff3\n";
		die "Error, cannot find $annot_gff3" unless (-s $annot_gff3);
		open my $fh, '<', $annot_gff3 or die "Error, cannot read $annot_gff3";
		GFF: while (my $line=<$fh>) {
        	chomp $line;

			# Remove carriage returns
			$line =~ s/\r//g;

			# Find features of interest
			next GFF if($line =~ m/^$/);
			next GFF if($line =~ m/^#/);
			my @cols = split "\t", $line;
			my ($source, $type, $gene_info) = ($cols[1], $cols[2], $cols[8]);
			next GFF if(!defined $type);
			next GFF unless ($type eq $feature);
			$cols[1] = $genome;

			# Find parent
			my @description_parts = split /$desc_seperator/, $gene_info;
			my $feature_parents = $description_parts[$desc_column];
			if($count eq 0) { warn "combine_all_gff3_files_in_repo: removing $desc_replace from IDs for $genome...\n"; }
			$feature_parents =~ s/$desc_replace//g;
			if($count eq 0) { warn "combine_all_gff3_files_in_repo: saving IDs (e.g. $feature_parents) for $genome...\n"; }

			# Check for duplicate features in same genome
			if(defined $ids_in_one_genome{$feature_parents}) {
				if($duplicate_id_count eq 0) { warn "combine_all_gff3_files_in_repo: found duplicate IDs in same genome. Excluding duplicates.\n"; }
				$duplicate_id_count = 1;
				next GFF;
			}
			$ids_in_one_genome{$feature_parents} = 1;

			# Check for unique transcript IDs
			if(defined $unique_transcript_ids{$feature_parents}) {

				# warn
				if($change_id_count eq 0) {
					warn "combine_all_gff3_files_in_repo: found IDs multiple times, so generating unique IDs\n";
					$change_id_count=1;
				}

				# Make a new ID
				my $new_transcript_id = $feature_parents . '_synima_' . time();
				while(defined $unique_transcript_ids{$new_transcript_id}) {
					$new_transcript_id = $feature_parents . '_synima_' . time();
				}
				$new_transcript_ids{$genome}{'old_to_new'}{$feature_parents} = $new_transcript_id;
				#$new_transcript_ids{$genome}{'new_to_old'}{$new_transcript_id} = $feature_parents;
				$feature_parents = $new_transcript_id;
			}
			$unique_transcript_ids{$feature_parents} = 1;

			# Print
			$cols[8] = $feature_parents;
			my $new_line = join "\t", @cols;
			print $ofh "$new_line\n";
			$count++;
		}
		close $fh;
		die "Error: Found no $feature features in $annot_gff3. Check and re-run\n" if($count eq 0);
		warn "combine_all_gff3_files_in_repo: found $count $feature features in $annot_gff3\n";
	}
	close $ofh;
	return \%new_transcript_ids;
}

sub save_genome_codes_from_synima_parsed_gff3 {
	my $annots_gff3 = $_[0];
	my (%trans_id_to_genome, %genome_to_code);

	# Open GFF and save ID->Genome
	my $genome_seen;
	my $count = 0;
	open my $fh, '<', $annots_gff3 or die "Cannot open $annots_gff3 : $!\n";
	warn "Saving ID->Genome from $annots_gff3...\n";
	GFF: while (my $line=<$fh>) {
		chomp $line;
        	next GFF unless ($line =~/\w/);
        	next GFF if($line =~ m/^#/);
        	my @x = split /\t/, $line;
		my ($genome, $type, $gene) = ($x[1], $x[2], $x[8]);
		die "Corrupt $annots_gff3: Re-run Create_full_repo_sequence_databases.pl\n" if(!defined $gene);
		if((!defined $genome_seen) || ($genome_seen ne $genome)) { 
			$count = 0; 
			$genome_seen = $genome;
			$genome_to_code{$genome} = 1; # just store for now.
		}
		$trans_id_to_genome{$gene} = $genome;
	}
	close $fh;
	return (\%trans_id_to_genome, \%genome_to_code);
}

sub save_gene_struct_from_gff {
	my ($gff_file, $organism) = @_;

	# Open GFF and save struct
	my %genes;
	open my $fh, '<', $gff_file or die "Error, cannot read file $gff_file : $!\n";
	warn "save_gene_struct_from_gff: parsing $gff_file for $organism\n";
	GFF: while (my $line=<$fh>) {
		chomp $line;
		next unless ($line =~/\w/);
		next if ($line =~/^\#/);

		my @x = split /\t/, $line;
		my ($contig, $source, $type, $lend, $rend, $score, $orient, $phase, $gene) = @x;
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		# Include org in the gene identifier name. Also have ID as alias now.
		my $id = $organism . ":" . $gene;
		my $struct = &gff_line_make_struct($organism, $contig, $lend, $rend, $end5, $end3, $orient, $id, $id, $source, $type);
		# Check if its already been seen
		#die "Error, already stored info for a feature with ID = $id - Check GFF and re-run. ($line)\n" if (exists $genes{$id});
		$genes{$id} = $struct;
	}
	return \%genes;
}

sub gff_line_make_struct {
	my ($organism, $contig, $lend, $rend, $end5, $end3, $orient, $id, $alias, $source, $type) = @_;
	my $struct = { 
		org => $organism,
		source => $source,
		contig => $contig,
		feature => $type,
		orgMol => join(";", ($organism, $contig)),
		lend => $lend,
		rend => $rend,
		mid => int(($lend + $rend)/2),
		end5 => $end5,
		end3 => $end3,
		pos => undef,  # assign later as relative position of gene on the contig.
		orient => $orient,
		gene_id => $id,
		alias => $alias,
		mRNAs => [],
	};
	return $struct;
}

sub split_gff_dictionary_by_species {
	my ($input_gff, $repo_spec) = @_;
	die "Cannot find $input_gff. Re-run Create_full_repo_sequence_databases.pl\n" if(! -e $input_gff);

	# Parse repo spec file
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();

	# Input seq dictionary
	my $gff_struct = &save_gene_struct_from_gff($input_gff, 'all');

	# Split GFF by genomes
	warn "Splitting GFF by genomes...\n";
	foreach my $genome (@genomes) {
		my $output_file = $data_manager->get_data_dump_filename($genome, "Annotation");
		$output_file .= ".synima-parsed.GFF3";
		warn "$genome = $output_file...\n";
		open my $ofh, '>', $output_file or die "Cannot open output file $output_file : $!\n";
		foreach my $transcript_id(keys %{$gff_struct}) {
			my $genome_found = $$gff_struct{$transcript_id}{'source'};
			next if($genome_found ne $genome);
			my $contig  = $$gff_struct{$transcript_id}{'contig'};
			my $feature = $$gff_struct{$transcript_id}{'feature'};
			my $start   = $$gff_struct{$transcript_id}{'lend'};
			my $stop    = $$gff_struct{$transcript_id}{'rend'};
			my $strand  = $$gff_struct{$transcript_id}{'orient'};
			my $gene_id = $$gff_struct{$transcript_id}{'gene_id'};
			#$gene_id =~ s/^all:/$genome:/;
			$gene_id =~ s/^all://;
			print $ofh "$contig\t$genome\t$feature\t$start\t$stop\t.\t$strand\t.\t$gene_id\n";
		}
		close $ofh;
	}
	return 1;
}

1;
