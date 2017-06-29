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
	my ($input, $feature_wanted, $desc_seperator, $desc_column, $desc_replace) = @_;
	my (%hash_info, %hash_strand);
	warn "gff_to_contig_parent_to_cds_hash: saving $feature_wanted from $input (split col $desc_column by $desc_seperator and remove $desc_replace)...\n";
	open IN, "<$input" or die "Cannot open $input: $!\n";
	GFF: while(my $line = <IN>) {
		chomp $line;
		my @bits = split /\t/, $line;
		next GFF if ((@bits < 8) || ($line =~ m/^\#/));
		my ($contig, $source, $feature, $start, $stop, $score, $strand, $frame, $full_description) = @bits;
		next GFF if(($feature ne $feature_wanted) && ($feature ne 'all'));
		
		# Find parent
		my @description_parts = split /$desc_seperator/, $full_description;
		my $feature_parents = $description_parts[$desc_column];
		$feature_parents =~ s/$desc_replace//g;
		my $CDS = "$start-$stop";

		# Save (avoid leaving a space at the end)
		if(!defined $hash_info{$contig}{$feature_parents}) { $hash_info{$contig}{$feature_parents} = $CDS; }
		else { $hash_info{$contig}{$feature_parents} .= " $CDS"; }
		$hash_strand{$feature_parents} = $strand;
	}
	close IN;
	return (\%hash_info, \%hash_strand);
}

sub combine_all_gff3_files_in_repo {
	my ($repo_spec, $output) = @_;
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();

	warn "Combine_all_gff3_files_in_repo: printing to $output\n";
	die "$output already exists. Delete then re-run\n" if(-e $output);
	open my $ofh, '>>', $output or die "Error, cannot open file $output\n";
	foreach my $genome (@genomes) {
        	my $annot_gff3 = $data_manager->get_data_dump_filename($genome, "Annotation");
        	die "Error, cannot find $annot_gff3" unless (-s $annot_gff3);

		open my $fh, '<', $annot_gff3 or die "Error, cannot read $annot_gff3";
		GFF: while (my $line=<$fh>) {
        		chomp $line;
			next GFF if($line =~ m/^$/);
        		my @cols = split "\t", $line;
        		$cols[1] = $genome;
        		my $new_line = join "\t", @cols;
            		print $ofh "$new_line\n" ;
		}
		close $fh
	}
	close $ofh;
	return 1;
}

sub save_genome_codes_from_gff3 {
	my ($annots_gff3, $feature, $desc_seperator, $desc_column, $desc_replace) = @_;
	my (%trans_id_to_genome, %genome_to_code);

	# Open GFF and save ID->Genome
	my $count = 0;
	open my $fh, '<', $annots_gff3 or die "Cannot open $annots_gff3 : $!\n";
	warn "Saving ID->Genome from $annots_gff3...\n";
	GFF: while (my $line=<$fh>) {
		chomp $line;
        	next GFF unless ($line =~/\w/);
        	my @x = split /\t/, $line;
		my ($genome, $type, $gene_info) = ($x[1], $x[2], $x[8]);
		next GFF unless ($type eq $feature);

		# Broad Style:
		#if(($gene_info =~ m/ID=([^;]+)/) && ($gene_info =~ /Alias=([^;]+)/) && ($gene_info =~ /Parent=([^;]+)/)) {
		if(($gene_info =~ m/ID=([^;]+)/) && ($gene_info =~ /Parent=([^;]+)/)) {
			# Find transcript ID
			my $trans_id;
        		if ($gene_info =~ /ID=([^;]+)/) {
        			$trans_id = $1;
				$trans_id_to_genome{$trans_id}= $genome;
			} else { die "Error, cannot find trans_id from $line\n"; }

			my $alias_id;
			if ($gene_info =~ /Alias=([^;]+)/) {
				$alias_id = $1;
				$trans_id_to_genome{$alias_id}= $genome;
			}
			if($count eq 0) { warn "Saving ID (e.g. $trans_id) and optionally Alias ($alias_id) from GFF (Broad format). If these ID's are not present in your fasta files - quit and change settings...\n"; $count = 1; }
		} 

		# 1 ID
		else {

			# Find parent
			my @description_parts = split /$desc_seperator/, $gene_info;
			my $feature_parents = $description_parts[$desc_column];
			if($count eq 0) { warn "trying to remove $desc_replace\n"; }
			$feature_parents =~ s/$desc_replace//g;
			if($count eq 0) { warn "Saving ID (e.g. $feature_parents) to genome (e.g. $genome) from GFF (Normal format). If these ID's are not present in your fasta files - quit and change settings...\n"; $count = 1; }
			$trans_id_to_genome{$feature_parents} = $genome;
		}
		$genome_to_code{$genome} = 1; # just store for now.
	}
	close $fh;
	return (\%trans_id_to_genome, \%genome_to_code);
}

sub save_gene_struct_from_gff {
	my ($gff_file, $organism, $feature, $desc_seperator, $desc_column, $desc_replace) = @_;

	my %genes;

	# Open GFF and save ID->Genome
	my $count = 0;
	open my $fh, '<', $gff_file or die "Error, cannot read file $gff_file : $!\n";
	warn "\t-parsing $gff_file for $organism\n";
	GFF: while (my $line=<$fh>) {
		chomp $line;
		next unless ($line =~/\w/);
		next if ($line =~/^\#/);

		my @x = split /\t/, $line;
		my ($contig, $source, $type, $lend, $rend, $score, $orient, $phase, $gene_info) = @x;
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		#next unless ($feat_type =~ /gene|mrna/i);
		next GFF unless ($type eq $feature);

		# Broad Style:
		if(($gene_info =~ m/ID=([^;]+)/) && ($gene_info =~ /Alias=([^;]+)/) && ($gene_info =~ /Parent=([^;]+)/)) {

			# Include org in the gene identifier name:
			$gene_info =~ /ID=([^;\s]+)/;
			my $id = $1 or die "read_GFF: Error, no ID from $gene_info - Check settings and re-run.\n";
			$id = ($organism . ":" . $id);
			my $alias;

			if ($gene_info =~ /Alias=([^\;\s]+)/) { $alias = $1; }

			if ($type eq 'gene') {
				my $struct = &gff_line_make_struct($organism, $contig, $lend, $rend, $end5, $end3, $orient, $id, $alias);

				# Check if its already been seen
				die "read_GFF: Error, already stored info for a feature with ID = $id" if (exists $genes{$id});
				$genes{$id} = $struct;
			}
			elsif ($type eq 'mRNA') {
				$gene_info =~ /Parent=([^\;\s]+)/ or die "read_GFF: Error, cannot extract parent from $gene_info - Check settings and re-run.\n";
				my $gene_id = $1 or die "read_GFF: Error, no gene ID extracted from mRNA $gene_info - Check settings and re-run.\n";
				$gene_id = $organism . ":" . $gene_id; # remember to include organism.

				my $gene_struct = $genes{$gene_id} or die "read_GFF: Error, no gene stored for $gene_id - Check settings and GFF and re-run.\n";
				push (@{$gene_struct->{mRNAs}}, $id); # store mRNAs linked to gene.
			} else {
				# Save something
				my $struct = &gff_line_make_struct($organism, $contig, $lend, $rend, $end5, $end3, $orient, $id, $alias);

				# Check if its already been seen
				die "Error, already stored info for a feature with ID = $id  - Check GFF and re-run.\n" if (exists $genes{$id});
				$genes{$id} = $struct;
			}
		}
		# 1 ID
		else {
			# Find ID
			my @description_parts = split /$desc_seperator/, $gene_info;
			my $feature_parents = $description_parts[$desc_column];
			if($count eq 0) { warn "-trying to remove $desc_replace\n"; }
			$feature_parents =~ s/$desc_replace//g;
			if($count eq 0) { warn "Saving ID (e.g. $feature_parents) to genome (e.g. $organism) from GFF (Normal format). If these ID's are not present in your fasta files - quit and change settings...\n"; $count = 1; }
			#$trans_id_to_genome{$feature_parents} = $genome;
			
			# Include org in the gene identifier name. Also have ID as alias now.
			my $id = $organism . ":" . $feature_parents;
			my $struct = &gff_line_make_struct($organism, $contig, $lend, $rend, $end5, $end3, $orient, $id, $id);

			# Check if its already been seen
			die "Error, already stored info for a feature with ID = $id - Check GFF and re-run.\n" if (exists $genes{$id});
			$genes{$id} = $struct;
		}
	}
	return \%genes;
}

sub gff_line_make_struct {
	my ($organism, $contig, $lend, $rend, $end5, $end3, $orient, $id, $alias) = @_;	
	my $struct = { 
		org => $organism,
		contig => $contig,
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

1;
