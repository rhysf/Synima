package daglines;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use Data::Dumper;

### rfarrer@broadinstitute.org

sub save_aligncoords {
	my $file = $_[0];
	my (%gene_synteny, %genes);
	$/ = "##";
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	BLOCK: while(my $alignment_block=<$fh>) {
		next BLOCK if ($alignment_block eq '##');
		my @lines = split /\n/, $alignment_block;

		# Alignment info
		GENES: for(my $i=1; $i<scalar(@lines); $i++) {
			next GENES if($lines[$i] =~ m/##/);
			my @alignment_info = split /\s+/, $lines[$i];
			my ($genome1, $chr1, $genome1_gene_name, $start1, $stop1, $ignore1, $ignore2, $genome2, $chr2, $genome2_gene_name, $start2, $stop2) = @alignment_info;
			die "Error: Unexpected aligncoords format for $file. Check settings/config file. line = $lines[$i]\n" if(!defined $stop2);

			# Sometimes genome IDs have a space followed by a description. Check there is a : in the expected variable. If not, remove first, and reinitalise array
			if(($genome1_gene_name !~ m/\:/) && ($start1 =~ m/\:/)) {
				splice @alignment_info, 2, 1;
				($genome1, $chr1, $genome1_gene_name, $start1, $stop1, $ignore1, $ignore2, $genome2, $chr2, $genome2_gene_name, $start2, $stop2) = @alignment_info;
			}
			if(($genome2_gene_name !~ m/\:/) && ($start2 =~ m/\:/)) {
				splice @alignment_info, 9, 1;
				($genome1, $chr1, $genome1_gene_name, $start1, $stop1, $ignore1, $ignore2, $genome2, $chr2, $genome2_gene_name, $start2, $stop2) = @alignment_info;
			}

			# Genomes and gene name
			my @gen1 = split /\:/, $genome1_gene_name;
			my @gen2 = split /\:/, $genome2_gene_name;
			die "Error: Unexpected aligncoords format for genome and gene name (1) separated by ; ($genome1_gene_name) on line $lines[$i]\n" if(!defined $gen1[1]);
			die "Error: Unexpected aligncoords format for genome and gene name (2) separated by ; ($genome2_gene_name) on line $lines[$i]\n" if(!defined $gen2[1]);

			# Save gene positions
			$genes{$genome1}{$chr1}{$gen1[1]}{$start1} = $stop1;
			$genes{$genome2}{$chr2}{$gen2[1]}{$start2} = $stop2;

			# Save gene synteny
			my $start_stop1 = $start1 . '-' . $stop1;
			my $start_stop2 = $start2 . '-' . $stop2;
			$gene_synteny{$genome1}{$genome2}{$chr1}{$chr2}{$start_stop1} = $start_stop2;
			$gene_synteny{$genome2}{$genome1}{$chr2}{$chr1}{$start_stop2} = $start_stop1;
		}
	}
	close $fh;
	$/ = "\n";
	die "Nothing saved from $file. Check file and that is confirms to save_aligncoords format.\n" if((! %gene_synteny) || (! %genes));	
	return (\%gene_synteny, \%genes);
}

sub save_aligncoords_spans {
	my ($file) = $_[0];
	my (%spans_info);
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	#warn "Saving synteny from $file...\n";
	while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t|\s+/, $line;
		my ($genome_and_contig1, $ss1, $size1, $genome_and_contig2, $ss2, $size2) = (@bits);
		die "unexpected aligncoords.spans format for $file. Check settings/config file. line = $line\n" if(!defined $ss2);
		my @gc1 = split /;/, $genome_and_contig1;	
		my @gc2 = split /;/, $genome_and_contig2;
		$spans_info{$gc1[0]}{$gc2[0]}{$gc1[1]}{$gc2[1]}{$ss1} = $ss2;
		$spans_info{$gc2[0]}{$gc1[0]}{$gc2[1]}{$gc1[1]}{$ss2} = $ss1;
	}
	close $fh;
	die "Nothing saved from $file. Check file and that is confirms to save_aligncoords.spans format.\n" if(! %spans_info);	
	warn "save_aligncoords_spans: saved info for " . scalar(keys(%spans_info)) . " genomes in $file\n";
	return (\%spans_info);
}

sub save_genome_names_hash_from_aligncoords_spans {
	my ($file) = $_[0];
	my (%spans_info);
	open my $fh, '<', $file or die "Cannot open $file: $!\n";
	warn "Saving genome names from $file...\n";
	SPANS: while(my $line=<$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($genome_and_contig1, $ss1, $size1, $genome_and_contig2, $ss2, $size2) = (@bits);
		die "Unpected format of $line in file $file\n" if(!defined $size2);

		my @gc1 = split /;/, $genome_and_contig1;
		my @gc2 = split /;/, $genome_and_contig2;
		die "Unexpected format of $genome_and_contig1 (expected ; between genome and contig) on line $line in file $file\n" if(!defined $gc1[1]);
		die "Unexpected format of $genome_and_contig2 (expected ; between genome and contig) on line $line in file $file\n" if(!defined $gc2[1]);
		$spans_info{$gc1[0]} = 1;
		$spans_info{$gc2[0]} = 1; 
	}
	close $fh;
	die "Nothing saved from $file. Check file and that is confirms to save_aligncoords.spans format.\n" if(! %spans_info);	
	return (\%spans_info);
}

sub save_contig_lengths_from_config {
	my $config_data = $_[0];
	my %contig_lengths;
	CONFIG: foreach my $config_keys(keys %{$config_data}) {
		next CONFIG unless($config_keys =~ m/^genome_contig_lengths_/);
		my $genome = $config_keys;
		$genome =~ s/^genome_contig_lengths_//;
		$contig_lengths{$genome} = $$config_data{$config_keys};
	}
	die "Nothing saved in save_contig_lengths_from_config.\n" if(! %contig_lengths);	
	return \%contig_lengths;
}

sub save_aligncoords_and_reverse_wrapper_from_config {
	my $config_data = $_[0];
	warn "save_aligncoords_and_reverse_wrapper_from_config...\n";
	my ($gene_synteny, $genes) = &save_aligncoords($$config_data{'aligncoords_location'}); 
	$gene_synteny = &reverse_placement_on_chromosome_for_genes_hash($gene_synteny, $config_data);
	warn "save_aligncoords_and_reverse_wrapper_from_config: finished\n";
	return ($gene_synteny, $genes);
}

sub save_aligncoords_spans_and_reverse_wrapper_from_config {
	my ($config_data, $file_name) = @_;
	warn "save_aligncoords_spans_and_reverse_wrapper_from_config...\n";
	my $genome_synteny = &save_aligncoords_spans($$config_data{$file_name});
	$genome_synteny = &reverse_placement_on_chromosome_for_genomes_hash($genome_synteny, $config_data);
	return $genome_synteny;
}

sub genes_subset_from_gene_list_from_config {
	my ($config_data, $file_name, $genes) = @_;

	# Save column 1
	$/ = "\n";
	my $gene_ids_want = tabfile::save_one_column($$config_data{$file_name}, 0);

	# Find gene IDs in genes hash and save
	my %genes_subset;
	foreach my $genome(keys %{$genes}) {
		foreach my $chr(keys %{$$genes{$genome}}) {
			GENES: foreach my $gene_ids(keys %{$$genes{$genome}{$chr}}) {
				next GENES unless(defined $$gene_ids_want{$gene_ids});
				foreach my $start(keys %{$$genes{$genome}{$chr}{$gene_ids}}) {
					my $stop = $$genes{$genome}{$chr}{$gene_ids}{$start};
					$genes_subset{$genome}{$chr}{$gene_ids}{$start} = $stop;
				}
			}
		}
	}
	return (\%genes_subset);
}

sub genome_to_chromosome_to_length_hash_subset {
	my $config_data = $_[0];
	#warn "genome_to_chromosome_to_length_hash_subset...\n";
	my (%chr_reversed);
	return (\%chr_reversed) if(!defined $$config_data{'reverse_compliment'});

	# Save contig lengths
	my $genome_contig_lengths = &save_contig_lengths_from_config($config_data);

	# Save genome->chromosome hash subset of those to be reversed
	my @chr_split = split /,/, $$config_data{'reverse_compliment'};
	foreach my $genome_and_contigs(@chr_split) { 
		my @genome_and_chr = split /:/, $genome_and_contigs;
		die "No genome $genome_and_chr[0] found in config (reverse_compliment). Check or remake\n" if(!defined $$genome_contig_lengths{$genome_and_chr[0]});
		die "No contig specified for genome $genome_and_chr[0] in config (reverse_compliment). Check or remake\n" if(!defined $genome_and_chr[1]);
		die "no length found for genome $genome_and_chr[0] contig $genome_and_chr[1] in config (reverse_compliment). Check or remake\n" if(!defined $$genome_contig_lengths{$genome_and_chr[0]}{$genome_and_chr[1]});
		$chr_reversed{$genome_and_chr[0]}{$genome_and_chr[1]} = $$genome_contig_lengths{$genome_and_chr[0]}{$genome_and_chr[1]};
	}
	return (\%chr_reversed);
}

sub check_hash_keys_from_aligncoords {
	my ($aligncoords_spans_hash, $number_of_keys, $key1, $key2, $key3, $key4, $key5) = @_;
	my $found = 0;

	# Check individiually, otherwise Perl creates keys to check for subkeys
	if(defined $$aligncoords_spans_hash{$key1}) {
		if(defined $$aligncoords_spans_hash{$key1}{$key2}) { 
			if(defined $$aligncoords_spans_hash{$key1}{$key2}{$key3}) {
				if(defined $$aligncoords_spans_hash{$key1}{$key2}{$key3}{$key4}) { 
					if($number_of_keys eq 4) { $found = 1; } 
					else {
						if(defined $$aligncoords_spans_hash{$key1}{$key2}{$key3}{$key4}{$key5}) { $found = 1; }
					}
				}
			}
		}
	}
	return $found;
}

sub reverse_placement_on_chromosome_for_genomes_hash {
	my ($ref_to_hash_to_reverse, $config_data) = @_;
	warn "reverse_placement_on_chromosome_for_genomes_hash...\n";

	# Save genome->chromosome hash subset of those to be reversed
	my $chr_reversed = &genome_to_chromosome_to_length_hash_subset($config_data);
	return $ref_to_hash_to_reverse if (! %{$chr_reversed});

	# Genes, wAnnots are like so: my $gene_stop = $$genes{$name1}{$chromosomes}{$gene_start};
	foreach my $genome(keys %{$ref_to_hash_to_reverse}) {
		#warn "$genome..\n";
		CHR: foreach my $chr(keys %{$$ref_to_hash_to_reverse{$genome}}) {
			#warn "$chr ?\n";
			next CHR if(!defined $$chr_reversed{$genome}{$chr});
			#warn "what is in my $genome $chr ?\n";
			foreach my $gene_start(keys %{$$ref_to_hash_to_reverse{$genome}{$chr}}) {
				my $gene_stop = $$ref_to_hash_to_reverse{$genome}{$chr}{$gene_start};

				# delete old entry
				delete $$ref_to_hash_to_reverse{$genome}{$chr}{$gene_start};

				# make new entry
				# stop and start are purposefully switched here
				my $new_start = $$chr_reversed{$genome}{$chr} - $gene_stop;
				my $new_stop = $$chr_reversed{$genome}{$chr} - $gene_start;
				$$ref_to_hash_to_reverse{$genome}{$chr}{$new_start} = $new_stop;
				#warn "$gene_start -> $new_start and $gene_stop -> $new_stop\n";
			}
		}
	}
	return $ref_to_hash_to_reverse;
}

sub reverse_placement_on_chromosome_for_genes_hash {
	my ($ref_to_hash_to_reverse, $config_data) = @_;

	# Save genome->chromosome hash subset of those to be reversed
	my $chr_reversed = &genome_to_chromosome_to_length_hash_subset($config_data);
	return $ref_to_hash_to_reverse if (! %{$chr_reversed});

	# Gene synteny is like so: my $gene_start_stop2(keys %{$$gene_synteny{$name1}{$name2}{$chromosomes}{$chr2}{$gene_start_stop1}}
	foreach my $genome(keys %{$ref_to_hash_to_reverse}) {
		GENOME2: foreach my $genome2(keys %{$$ref_to_hash_to_reverse{$genome}}) {
			next GENOME2 if((!defined $$chr_reversed{$genome}) && (!defined $$chr_reversed{$genome2}));
			foreach my $chr1(keys %{$$ref_to_hash_to_reverse{$genome}{$genome2}}) {
				CHRS: foreach my $chr2(keys %{$$ref_to_hash_to_reverse{$genome}{$genome2}{$chr1}}) {
					# neither want flipping
					next CHRS if((!defined $$chr_reversed{$genome}{$chr1}) && (!defined $$chr_reversed{$genome2}{$chr2}));

					# At least 1 of them needs flipping
					foreach my $chr1_ss(keys %{$$ref_to_hash_to_reverse{$genome}{$genome2}{$chr1}{$chr2}}) {
						my $chr2_ss = $$ref_to_hash_to_reverse{$genome}{$genome2}{$chr1}{$chr2}{$chr1_ss};
						die "start and stop not defined: 1) $chr1_ss and 2) $chr2_ss\n" if((!defined $chr1_ss) || (!defined $chr2_ss));

						my @chr1ss = split /-/, $chr1_ss;
						my @chr2ss = split /-/, $chr2_ss;
						foreach(@chr1ss, @chr2ss) { die "start or stop not defined in genes hash ($chr1_ss -> $chr2_ss): $_\n" if(!defined $_); } 
						delete $$ref_to_hash_to_reverse{$genome}{$genome2}{$chr1}{$chr2}{$chr1_ss};

						# just genome 1 chr. needs flipping
						if((defined $$chr_reversed{$genome}{$chr1}) && (!defined $$chr_reversed{$genome2}{$chr2})) {
							my $new_start1 = $$chr_reversed{$genome}{$chr1} - $chr1ss[1];
							my $new_stop1 = $$chr_reversed{$genome}{$chr1} - $chr1ss[0];
							my $new_ss1 = ($new_start1 . '-' . $new_stop1);
							$$ref_to_hash_to_reverse{$genome}{$genome2}{$chr1}{$chr2}{$new_ss1} = $chr2_ss;
						}
						# just genome 2 chr. needs flipping
						elsif((! defined $$chr_reversed{$genome}{$chr1}) && (defined $$chr_reversed{$genome2}{$chr2})) {
							my $new_start2 = $$chr_reversed{$genome2}{$chr2} - $chr2ss[1];
							my $new_stop2 = $$chr_reversed{$genome2}{$chr2} - $chr2ss[0];
							my $new_ss2 = ($new_start2 . '-' . $new_stop2);
							$$ref_to_hash_to_reverse{$genome}{$genome2}{$chr1}{$chr2}{$chr1_ss} = $new_ss2;	
						}
						# both
						else {
							if(!defined $$chr_reversed{$genome}{$chr1}) { warn "chr_reversed{$genome2}{$chr1} not defined (should be both defined?)\n"; }
							if(!defined $$chr_reversed{$genome2}{$chr2}) { warn "chr_reversed{$genome2}{$chr2} not defined (should be both defined?)\n"; }
							my $new_start1 = $$chr_reversed{$genome}{$chr1} - $chr1ss[1];
							my $new_stop1 = $$chr_reversed{$genome}{$chr1} - $chr1ss[0];
							my $new_ss1 = ($new_start1 . '-' . $new_stop1);
							my $new_start2 = $$chr_reversed{$genome2}{$chr2} - $chr2ss[1];
							my $new_stop2 = $$chr_reversed{$genome2}{$chr2} - $chr2ss[0];
							my $new_ss2 = ($new_start2 . '-' . $new_stop2);
							$$ref_to_hash_to_reverse{$genome}{$genome2}{$chr1}{$chr2}{$new_ss1} = $new_ss2; 
						}
					}
				}
			}
		}
	}
	return $ref_to_hash_to_reverse;
}

sub find_contig_order_from_aligncoords_spans_and_fastas {
	my ($fasta_order1, $fasta_order2, $isolate_name1, $isolate_name2, $aligncoords_spans, $verbose) = @_;
	my @new_order2;
	my %contigs_found;
	ORDER: foreach my $chr_order(@{$fasta_order1}) {
		my $search_for = "grep \$'$isolate_name1;$chr_order\\t' $aligncoords_spans | grep \$'$isolate_name2;'";
		if($verbose eq 'y') { warn "Looking for syntenic regions between $isolate_name1 ($chr_order) and $isolate_name2... using $search_for\n"; }
		my $lines = `$search_for`;
		die "found nothing for: grep \$'$isolate_name1;$chr_order\t' $aligncoords_spans | grep $isolate_name2.... investigate.\n" if(!defined $lines);
		#next ORDER if(!defined $lines);
		my @line = split /\n/, $lines;
		my (%new_order_info, %align_lengths);
		foreach my $chr_line(@line) {
			my @bits = split /\t/, $chr_line;
			my ($is1, $ss1, $length1, $is2, $ss2, $length2) = @bits;
			my @isolate_and_chr1 = split /;/, $is1;
			my @isolate_and_chr2 = split /;/, $is2;
			my @start_stop1 = split /-/, $ss1;
			my @start_stop2 = split /-/, $ss2;

			# find order
			# I want order from isolate1, and then to save isolate2 chromosome names. 
			if($isolate_and_chr1[0] eq $isolate_name1) { $new_order_info{$start_stop1[0]}{$start_stop1[1]} = $isolate_and_chr2[1]; }
			if($isolate_and_chr2[0] eq $isolate_name1) { $new_order_info{$start_stop2[0]}{$start_stop2[1]} = $isolate_and_chr1[1]; }

			# Length (can be accumulative here) should be used to decide if multiple hits are found
			$align_lengths{$isolate_and_chr1[1]} += $length1;
			$align_lengths{$isolate_and_chr2[1]} += $length1;
		}

		# sort them and save order. Still has redundancies
		foreach my $start(sort { $a <=> $b } keys %new_order_info) {
			foreach my $stop(keys %{$new_order_info{$start}}) {
				my $contig = $new_order_info{$start}{$stop};

				# check if it is already there?
				if(defined $contigs_found{$contig}) {

					# if the last time it was found, it was smaller - then replace with this one
					if($contigs_found{$contig} < $align_lengths{$contig}) {
						# delete it from array
						for(my $i=0; $i<scalar(@new_order2); $i++) {
							if($new_order2[$i] eq $contig) { splice @new_order2, $i, 1; }
						}
						#warn "$contig moved to front due to $contigs_found{$contig} < $align_lengths{$contig}\n";

						# replace with new value
						push @new_order2, $contig;
						$contigs_found{$contig} = $align_lengths{$contig};
					}
				} else {
					my $length = $align_lengths{$contig};
					$contigs_found{$contig} = $align_lengths{$contig};
					push @new_order2, $contig;
					#warn "putative order (based on $chr_order): $contig\n";
				}
			}
		}
	}

	# Add missing chromosomes/contigs etc.
	my %found;
	my $new_order = "contig_order_$isolate_name2=";
	foreach(@new_order2) { 
		$new_order .= "$_,";
		$found{$_}=1; 
	}
	ORDER: foreach my $chr_order(@{$fasta_order2}) {
		if(!defined $found{$chr_order}) { $new_order .= "$chr_order,"; }
	}
	$new_order =~ s/,$//;

	# Return string for config, and new array
	return ($new_order, \@new_order2);
}

1;
