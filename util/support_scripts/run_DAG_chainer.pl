#!/usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Std;
use Storable qw(dclone);
use FindBin qw($Bin);
use lib "$Bin/../../modules/";
use read_GFF;
use Carp;
#use lib ("$FindBin::Bin/PerlLib");
use IniReader;
$|=1;

# Opening commands
my $usage = "Usage: perl $0 -c <Config file>
Optional: -z Only write DAGchainer format for inputs (y/n) [n]
	  -f Feature wanted from GFF [mRNA]
	  -s Seperator in GFF description for gene names (\" ; etc) [;]
	  -d GFF description part number with the parent/gene info [0]
	  -m Remove additional comments in column [Parent=]
	  -v Verbose (y/n) [n]
Notes:    FASTA ID's and GFF ID's must match up 
          Use settings -f, -s, -d and -m to ensure they do\n";
our($opt_c, $opt_d, $opt_f, $opt_m, $opt_s, $opt_v, $opt_z);
getopt('cdfmsvz');
die $usage unless($opt_c);
if(!defined $opt_d) { $opt_d = 0; }
if(!defined $opt_f) { $opt_f = 'mRNA'; }
if(!defined $opt_m) { $opt_m = 'Parent='; }
if(!defined $opt_s) { $opt_s = ';'; }
if(!defined $opt_v) { $opt_v = 'n'; }
if(!defined $opt_z) { $opt_z = 'n'; }
die "-d needs to be a number: $opt_d\n" unless($opt_d =~ m/^\d+$/);

# Init
my $DEBUG = 0;
my $uname = $ENV{HOSTTYPE};
my $progpath = "$FindBin::Bin/";
warn "Program path: $0\n";
warn "Settings: -c $opt_c -z $opt_z -f $opt_f -s $opt_s -d $opt_d -m $opt_m -v $opt_v\n";

# Save dagchainer config file
my ($dag_config, $organism_to_gff_file) = &save_dagchainer_config($opt_c);

# hidden opt for troubleshooting.
our $KEEP_DELCHER_FILES = 0; 

main: {
	 
	## Parse input match data, populate data structures:
	$|++;

	my $orthologs;
	if ($$dag_config{'orthologs_file'}) { $orthologs = &parse_orthologs($$dag_config{'orthologs_file'}); }
	
	if($opt_v eq 'y') { warn "-parsing matches file: $$dag_config{'inputFile'}\n"; }
	my @match_pairs = &parse_matches_file($$dag_config{'inputFile'}, $orthologs, $$dag_config{'max_e_value'});
	
	if($opt_v eq 'y') { warn "-parsing gene features...\n"; }
	my @gene_features = &parse_gff_files($organism_to_gff_file, $opt_f, $opt_s, $opt_d, $opt_m);

	# Assign lookup of feature based on gene_id, mRNA_id, and alias if available.
	if($opt_v eq 'y') { warn "-linking feature IDs to genes\n"; }
	my $feature_ID_to_gene = &link_feature_IDs_to_gene(@gene_features);
	
	# Order genes on scaffold and set relative position values.
	if($opt_v eq 'y') { warn "-ordering genes on scaffolds\n"; }
	my %orgMol_to_geneList = &splay_genes_on_molecule(@gene_features);

	# Get matches of gene pairs according to molecule comparisons.
	# Molecules and matches are sorted lexically
	# OrgMolPair is a tab-delimited key of orgMolA(tab)orgMolB
	if($opt_v eq 'y') { warn "-grouping matches by organism/molecule pairs\n"; }
	my %orgMolPair_to_matches = &group_matches_by_orgMolPairs(\@match_pairs, $feature_ID_to_gene);

	if ($$dag_config{'NOISE_FILTER_DIST'}) {
		if($opt_v eq 'y') { warn "-applying noise filter\n"; }
		%orgMolPair_to_matches = &apply_noise_filter(\%orgMolPair_to_matches, $feature_ID_to_gene, $$dag_config{'NOISE_FILTER_DIST'}, $$dag_config{'operation_mode'});
	}

	# open output file and force buffer flushing.
	my $aligncoords_out = ($$dag_config{'inputFile'} . '.aligncoords');
	if($opt_v eq 'y') { warn "-running DAGchainer -> $aligncoords_out\n"; }
	open my $ofh, '>', $aligncoords_out or die "Cannot open $aligncoords_out : $!\n"; 
	my $ref = select $ofh;
	$|++;
	select $ref;
	if($opt_v eq 'y') { print "Done parsing inputfile\n"; }

	# Perform the DAG chaining for each molecule pair (do +/- orientations separately)
	foreach my $molpair (sort keys %orgMolPair_to_matches) {
		print "molpair: $molpair\n" if (($DEBUG) || ($$dag_config{'SEE'}));
		my ($mol_1, $mol_2) = split /\t/, $molpair;

		next if (($mol_1 eq $mol_2) && (!$$dag_config{'INCLUDE_SELF'}));

		my $match_list_aref = $orgMolPair_to_matches{$molpair};

		if ($opt_z eq 'y') {
			&reformat_match_pair_dagchainer_output($mol_1, $mol_2, $match_list_aref, $feature_ID_to_gene, $ofh, $opt_v);
			next;
		}
		
		if($opt_v eq 'y') { print "***** Comparing $mol_1 to $mol_2 *****\n"; }
		## Create input file to cpp program.
		my $filename = "$$.delcher.input";
				
		my @pairIndexToAccs;
		open ART, ">$filename" or die "Can't open file.\n";
		my $pairID = 0;
		foreach my $match (@$match_list_aref) {
			#print "Match: " . Dumper($match);
			my ($acc_A, $acc_B, $e_value) = ($match->{accA}, $match->{accB}, $match->{e_value});
			if ($acc_A eq $acc_B) { die "Error, $acc_A eq B" . Dumper($match); }

			my $feature_A = $$feature_ID_to_gene{$acc_A} or die "Error, no feature retrieved based on $acc_A";
			my $feature_B = $$feature_ID_to_gene{$acc_B} or die "Error, no feature retreived based on $acc_B";

			my $midpt_A = $feature_A->{mid};
			my $midpt_B = $feature_B->{mid};
			if ($$dag_config{'operation_mode'} =~ /relative_position/i) {
				$midpt_A = $feature_A->{pos};
				$midpt_B = $feature_B->{pos};
			}
			
			my $e_value = $match->{e_value};
			my $score = scoringF($e_value, $$dag_config{'MAX_MATCH_SCORE'}, $$dag_config{'CONSTANT_MATCH_SCORE'});
			
			printf ART ("%d\t%d\t%d\t%f\n", $pairID, $midpt_A, $midpt_B, $score);
			
			$pairIndexToAccs[$pairID] = [$feature_A->{gene_id}, $feature_B->{gene_id}, $e_value];
			$pairID++;
		}
		close ART;
		
		#print Dumper(\@pairIndexToAccs);
		
		## forward direction
		&run_DAG_chainer($feature_ID_to_gene, $mol_1, $mol_2, $filename, \@pairIndexToAccs, "", $dag_config, $ofh, $opt_v); 

		## revcomp mol2
		&run_DAG_chainer($feature_ID_to_gene, $mol_1, $mol_2, $filename, \@pairIndexToAccs, "-r", $dag_config, $ofh, $opt_v);
		
		unlink ($filename) unless ($KEEP_DELCHER_FILES); #remove tempfile.
		
	} #end of molecule.

	close $ofh;
	exit(0);
}

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}

sub print_alignment {
	my ($mol1, $mol2, $match_header, $align_list_aref, $feature_ID_to_gene_href, $ofh, $verbose) = @_;
	my $ignore_alignment = 0;
	my $ignore_reason = "";

	my $num_aligned_pairs = scalar (@$align_list_aref);
	if($verbose eq 'y') { print "# $mol1 vs. $mol2 $match_header $num_aligned_pairs aligned pairs.\n"; }

	# more than one pair aligned.
	if (($opt_z eq 'n') && ($num_aligned_pairs < $$dag_config{'MIN_NUM_ALIGNED_PAIRS'})) { 
		$ignore_alignment = 1;
		$ignore_reason = " ignoring alignment.  Only $num_aligned_pairs, threshold=$$dag_config{'MIN_NUM_ALIGNED_PAIRS'}";
	}

	my $alignment_text = "## alignment $mol1 vs. $mol2 $match_header (num aligned pairs: $num_aligned_pairs):\n";
	my $IS_TANDEM = 0;
	my $IS_OVERLAPPING = 0;
	if (($opt_z eq 'n') && (!$ignore_alignment) && ($mol1 eq $mol2)) {
		my @acc1_coords;
		my @acc2_coords;
		print "## checking to see if should ignore a tandem aligment.\n" if ($$dag_config{'SEE'});
		## walk thru and see if any accession is showing up multiple times in alignment:
		my %seen;
		foreach my $alignedPair (@$align_list_aref) {
			my ($geneA_acc, $geneB_acc) = ($alignedPair->{acc_A}, $alignedPair->{acc_B});
			print "$geneA_acc,$geneB_acc\n" if($$dag_config{'SEE'});
			if (($geneA_acc ne '-' && $seen{$geneA_acc}) || ($geneB_acc ne '-' && $seen{$geneB_acc})) {
				print "SEEN!!, ignoring.\n" if ($$dag_config{'SEE'});
				$IS_TANDEM = 1;
				last;
			} else {
				$seen{$geneA_acc} = 1; 
				$seen{$geneB_acc} = 1;
			}
			## get coords for is_overlap determination.
			my ($feature_A, $feature_B) = ($feature_ID_to_gene_href->{$geneA_acc}, $feature_ID_to_gene_href->{$geneB_acc});
			
			my ($a_end5, $a_end3, $b_end5, $b_end3) = ($feature_A->{end5}, $feature_A->{end3}, $feature_B->{end5}, $feature_B->{end3});
			push (@acc1_coords, $a_end5, $a_end3);
			push (@acc2_coords, $b_end5, $b_end3);
		}

		@acc1_coords = sort {$a<=>$b} @acc1_coords;
		my $lend_A = shift @acc1_coords;
		my $rend_A = pop @acc1_coords;
		
		@acc2_coords = sort {$a<=>$b} @acc2_coords;
		my $lend_B = shift @acc2_coords;
		my $rend_B = pop @acc2_coords;

		#overlap
		if ($lend_A <= $rend_B && $rend_A >= $lend_B) { $IS_OVERLAPPING = 1; }
	}
    
	if (($opt_z eq 'n') && ($IS_TANDEM || $IS_OVERLAPPING)) {
		$ignore_alignment = 1;
		$ignore_reason = "ignoring tandem or overlapping self alignments.";
	} 

	if ($ignore_alignment && ($verbose eq 'y')) { print "# $ignore_reason\n"; } 
	else {
		print $ofh $alignment_text;
		print $alignment_text if ($$dag_config{'SEE'});
		
		foreach my $alignedPair (@$align_list_aref) {
			my ($acc_A, $acc_B, $dag_position_score, $Evalue) = ($alignedPair->{acc_A}, $alignedPair->{acc_B}, $alignedPair->{dag_position_score}, $alignedPair->{Evalue});
			my $feature_A = $feature_ID_to_gene_href->{$acc_A};
			my $feature_B = $feature_ID_to_gene_href->{$acc_B};
			my ($org_A, $mol_A, $id_A, $end5_A, $end3_A, $pos_A) = ($feature_A->{org}, $feature_A->{contig}, $feature_A->{acc}, $feature_A->{end5}, $feature_A->{end3}, $feature_A->{pos});
			my ($org_B, $mol_B, $id_B, $end5_B, $end3_B, $pos_B) = ($feature_B->{org}, $feature_B->{contig}, $feature_B->{acc}, $feature_B->{end5}, $feature_B->{end3}, $feature_B->{pos});
			my $outline = join("\t", $org_A, $mol_A, $acc_A,$end5_A, $end3_A, $pos_A, "MATCHES", $org_B, $mol_B, $acc_B, $end5_B, $end3_B, $pos_B, $Evalue, $dag_position_score) . "\n";
			print $ofh $outline;
			print $outline if ($$dag_config{'SEE'});
		}
	}
}

sub count_num_pairs {
    my $aref = shift;
    my $count = 0;
    foreach my $pair (@$aref) {
		my ($a, $b) = @$pair;
		#print "a: $a\tb: $b\n";
		unless ($a eq "-" || $b eq "-") {
			$count++;
		}
    }
    return ($count);
}

sub parse_matches_file {
	my ($file, $orthologs_href, $max_e_value) = @_;

	my %match_pairs;
	my %seenAcc; #track accessions examined already.
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while (my $line=<$fh>) {
		chomp $line;
		next unless ($line =~/\w/);
		my ($accA, $accB, $e_value, @rest) = split /\t/, $line;

		my $are_orthologs = &are_orthologous($accA, $accB, $orthologs_href);
		next if ((keys %{$orthologs_href}) && ($are_orthologs eq 0));

		# can't take logs of 0
		if ($e_value < 1e-250) { $e_value = 1e-250; }
		
		# run a few checks, see if we want this record.
		next if ($accA eq $accB);
		next unless ($e_value <= $max_e_value);
		
		($accA, $accB) = sort ($accA, $accB);			
		my $accPairKey =  join(";", ($accA, $accB));
		my $match = $match_pairs{$accPairKey};
		
		if (ref $match) {
			# Take the lowest e_value for the match pair.  Important when there are multiple HSPs reported between two accessions.
			if ($match->{e_value} > $e_value) { $match->{e_value} = $e_value; }
		} else {
			$match = { 
				accA => $accA,
				accB => $accB,
				e_value => $e_value 
			};
			$match_pairs{$accPairKey} = $match;
			
		}
	}
	close $fh;
	my @match_structs = values %match_pairs;
	return(@match_structs);
}

sub run_DAG_chainer {
	my ($feature_ID_to_gene_href, $mol1, $mol2, $filename, $pairIndexToAccs_aref, $reverseOrientFlag, $dag_config, $ofh, $verbose) = @_;

	# Run cpp program.
	if ($$dag_config{'SEE'}) {
		print "\n\n\nDAGCHAINERcpp input:\n";
		system "cat $filename";
	}

	my $tmpFile = ".$$.tmpOut";
	#my $cmd = "${progpath}dagchainer.$uname -G $GAP_LENGTH -O $GAP_OPEN_PENALTY -E $GAP_EXTENSION_PENALTY -S $MIN_ALIGNMENT_SCORE -D $MAX_DIST_BETWEEN_MATCHES  -F $filename $reverseOrientFlag > $tmpFile";
	my $cmd = "${progpath}dagchainer.x86_64-linux -G $$dag_config{'GAP_LENGTH'} -O $$dag_config{'GAP_OPEN_PENALTY'} -E $$dag_config{'GAP_EXTENSION_PENALTY'} -S $$dag_config{'MIN_ALIGNMENT_SCORE'} -D $$dag_config{'MAX_DIST_BETWEEN_MATCHES'} -F $filename $reverseOrientFlag > $tmpFile";
	#warn "run_DAG_chainer: $cmd \n";
	#print "CMD: $cmd \n";
	my $ret = system $cmd;
	die "ERROR, couldn't run command\n\n$cmd\n(ret: $ret)\n\nDid you recompile ${progpath}dagchainer.$uname for your OS?\n\n" if($ret);

	# Save DAGchainer output to memory
	my $all_alignments;
	open (OUT, "<$tmpFile");
	while (<OUT>) {
		print if($$dag_config{'SEE'});
		$all_alignments .= $_;
	}
	close OUT;
	unlink $tmpFile;

	# Print alignment
	my @individual_alignments = split (/>/, $all_alignments);
	shift @individual_alignments; #rid header that lacks alignment
	foreach my $individual_alignment (@individual_alignments) {
		my @align;
		my @matches = (split (/\n/, $individual_alignment));
		
		my $match_header = shift @matches;
		
		foreach my $match (@matches) {
			$match =~ s/^\s+//;
			
			if ($match =~ /^\d+:/) {
				chomp $match;
				my @x = split (/\s+/, $match);
				my ($index, $pairID, $pos1, $pos2, $match_score, $dag_chain_score) = @x;
				my $acc_pair_aref = $pairIndexToAccs_aref->[$pairID];
				my ($accA, $accB, $Evalue) = @$acc_pair_aref;
				push (@align, {acc_A=> $accA, acc_B=> $accB, dag_position_score => $dag_chain_score, Evalue => $Evalue, } );
			}
		}

		if ($reverseOrientFlag) { $match_header = "(reverse) $match_header"; }

		#print Dumper(\@align);
		&print_alignment($mol1, $mol2, $match_header, \@align, $feature_ID_to_gene_href, $ofh, $verbose);
	}
}

sub reformat_match_pair_dagchainer_output {
	my ($mol1, $mol2, $match_list_aref, $feature_ID_to_gene_href, $ofh, $verbose) = @_;
	my @dag_matches;
	foreach my $match_pair (@$match_list_aref) {
		my $struct = { 'acc_A'  => $match_pair->{accA},
				'acc_B'   => $match_pair->{accB},
				'Evalue'  => $match_pair->{e_value},
				'dag_position_score' => 0,
		};
		push (@dag_matches, $struct);
	}
	&print_alignment($mol1, $mol2, "", \@dag_matches, $feature_ID_to_gene_href, $ofh, $verbose);
}

sub scoringF {
	my ($evalue, $max_match, $constant_match) = @_;
	if ($constant_match) { return ($constant_match); } 
	
	my $matchScore = -log10($evalue);
	$matchScore *=10;
	$matchScore = int($matchScore +.5);
	$matchScore /= 10;
	$matchScore = $max_match if($matchScore > $max_match);
	return $matchScore;
}

sub parse_gff_files {
	my ($org_to_gff, $f, $s, $d, $m) = @_;
	
	my %genes;
	# assume read gene line before transcript line
	foreach my $organism (keys %{$org_to_gff}) {
		my $gff_file = $$org_to_gff{$organism};
		my $gene_for_org = gfffile::save_gene_struct_from_gff($gff_file, $organism, $f, $s, $d, $m);
		%genes = (%genes, %{$gene_for_org});
		#die "have i saved the right stuff?";
	}
	my @gene_list = values %genes;
	return(@gene_list);
}

sub link_feature_IDs_to_gene {
	my @genes = @_;

	my %features;
	foreach my $gene (@genes) {
		my $gene_id = $gene->{gene_id};
		my $alias = $gene->{alias};
		my @mRNA_ids = @{$gene->{mRNAs}};

		$features{$gene_id} = $gene;
		
		if ($alias) {
			die "Error, already stored feature based on alias: $alias" if ((exists $features{$alias}) && ($features{$alias} ne $gene));
			$features{$alias} = $gene;
		}

		foreach my $mRNA_id (@mRNA_ids) {
			die "Error, already stored a feature based on mRNA id: $mRNA_id" if (exists $features{$mRNA_id});
			$features{$mRNA_id} = $gene;
		}
	}
	return \%features;
}

sub splay_genes_on_molecule {
	my (@gene_features) = @_;

	my %orgMol_to_genes;
	foreach my $gene (@gene_features) {
		my $org = $gene->{org};
		my $mol = $gene->{contig};
		my $orgmol = $gene->{orgMol};
		push (@{$orgMol_to_genes{$orgmol}}, $gene);
	}

	# Assign relative position of genes
	foreach my $gene_list_aref (values %orgMol_to_genes) {
		my @genes = sort {$a->{mid}<=>$b->{mid}} @$gene_list_aref;
		my $counter = 0;
		foreach my $gene (@genes) {
			$counter++;
			$gene->{pos} = $counter;
		}
	}
	return(%orgMol_to_genes);
}

sub group_matches_by_orgMolPairs {
	my ($match_pairs_aref, $feature_ID_to_gene_href) = @_;

	my %orgMolPair_to_matches;
	foreach my $match_pair (@$match_pairs_aref) {
		my ($accA, $accB) = ($match_pair->{accA}, $match_pair->{accB});

		# no self matches!
		next if ($accA eq $accB);

		my $gene_struct_A = $feature_ID_to_gene_href->{$accA};
		my $gene_struct_B = $feature_ID_to_gene_href->{$accB};
		confess "Warning, cannot find a gene feature for [$accA]\n" unless ($gene_struct_A);
		confess "Warning, cannot find a gene feature for [$accB]\n" unless ($gene_struct_B);

		# Matches between two different isoforms of the same gene?  isoform self hits were already excluded on parsing.
		next if ($gene_struct_A eq $gene_struct_B);

		my $orgMolA = $gene_struct_A->{orgMol};
		my $orgMolB = $gene_struct_B->{orgMol};
		
		my $stored_match_pair = dclone($match_pair);

		# Keep lexically ordered
		if ($orgMolB lt $orgMolA) {
			# swap info:
			$stored_match_pair->{accA} = $match_pair->{accB};
			$stored_match_pair->{accB} = $match_pair->{accA},
		}
				
		my $orgMolPairKey = join ("\t", sort ($orgMolA, $orgMolB));
		
		push (@{$orgMolPair_to_matches{$orgMolPairKey}}, $stored_match_pair);
	}
	return(%orgMolPair_to_matches);
}

sub parse_orthologs {
	my ($orthologs_file) = @_;
	
	my %orthologs;
	open my $fh, '<', $orthologs_file or die "Error, cannot open $orthologs_file : $!\n";
	warn "Saving orthologs from $orthologs_file...\n";
	while (my $line=<$fh>) {
		chomp $line;
		my ($tag, $feats) = split /\t/, $line;
		$feats =~ s/^\s+|\s+$//g;
		my @genes = split(/\s+/, $feats);
		
		# strip parents off
		my @entries;
		foreach my $gene (@genes) {
			$gene =~ s/\(.*$//g;
			push (@entries, $gene);
		}

		# link to list of orthologs
		foreach my $entry (@entries) { $orthologs{$entry} = \@entries; }
	}
	close $fh;
	return (\%orthologs);
}

sub are_orthologous {
	my ($accA, $accB, $orthologs_href) = @_;

	if (my $list_aref = $orthologs_href->{$accA}) {
		
		foreach my $acc (@$list_aref) {
			if ($acc eq $accB) {
				return(1); # yes, orthologs
			}
		}
	}

	return(0);  # not orthologous
}

sub apply_noise_filter {
	my ($orgMolPair_to_matches_href, $feature_ID_to_gene_href, $noise_filt_dist, $operation_mode) = @_;

	my %refined_orgMolPair_to_matches;
	
	foreach my $orgMolPair (keys %$orgMolPair_to_matches_href) {
		
		my @matches = @{$orgMolPair_to_matches_href->{$orgMolPair}};

		my %geneA_to_matches;
		my %geneB_to_matches;
		foreach my $match (@matches) {
			my $accA = $match->{accA};
			my $accB = $match->{accB};
		   
			push (@{$geneA_to_matches{$accA}}, $match);
			push (@{$geneB_to_matches{$accB}}, $match);
		}

		foreach my $acc_A (keys %geneA_to_matches) {
			my @matches = @{$geneA_to_matches{$acc_A}};
			&filter_matches(\@matches, 'accB', $feature_ID_to_gene_href, $noise_filt_dist, $operation_mode);
		}

		foreach my $acc_B (keys %geneB_to_matches) {
			my @matches = @{$geneB_to_matches{$acc_B}};
			&filter_matches(\@matches, 'accA', $feature_ID_to_gene_href, $noise_filt_dist, $operation_mode);
		}

		foreach my $match (@matches) {
			if ($match->{best_in_filter}) {
				push (@{$refined_orgMolPair_to_matches{$orgMolPair}}, $match);
			}
		}
	}
	return(%refined_orgMolPair_to_matches);
}

sub filter_matches {
	my ($matches_aref, $gene_key, $feature_ID_to_gene_href, $noise_filt_dist, $operation_mode)  = @_;
	my @matches = @$matches_aref;
	@matches = sort { $feature_ID_to_gene_href->{ $a->{$gene_key} }->{pos} <=> $feature_ID_to_gene_href->{ $b->{$gene_key} }->{pos} } @matches;

	my $match = shift @matches;
	my @clusters = ([$match]);

	foreach my $match (@matches) {
		my $prev_cluster = $clusters[$#clusters];
		my $prev_match = $prev_cluster->[$#$prev_cluster];
		my $prev_gene_id = $prev_match->{$gene_key};
		my $prev_gene_struct =  $feature_ID_to_gene_href->{$prev_gene_id};
		confess "Error, no gene struct for $prev_gene_id" unless (ref $prev_gene_struct); 

		my $prev_coord = ($operation_mode =~ /position/i) ? $prev_gene_struct->{pos} : $prev_gene_struct->{mid};

		my $curr_gene_id =  $match->{$gene_key};
		my $curr_gene_struct = $feature_ID_to_gene_href->{$curr_gene_id};
		confess "Error, no gene struct for $curr_gene_id" unless (ref $curr_gene_struct);

		my $curr_coord = ($operation_mode =~ /position/i) ? $curr_gene_struct->{pos} :  $curr_gene_struct->{mid};

		# add to current cluster
		my $delta = ($curr_coord - $prev_coord);
		if ($delta <= $noise_filt_dist) { push (@$prev_cluster, $match); }

		# create a new cluster
		else { push (@clusters, [$match]); }
	}

	# Mark the best one in each cluster.
	foreach my $cluster (@clusters) {
		my @matches = @$cluster;
		@matches = sort {$a->{Evalue}<=>$b->{Evalue}} @matches;
		$matches[0]->{best_in_filter} = 1; # entry with lowest Evalue is the best match.
	}
	return;
}

sub save_dagchainer_config {
	my $conf = $_[0];

	my %config;
	my $properties_obj = new IniReader($conf);

	$config{'SEE'} = $properties_obj->get_value("Parameters", "Verbose");
	if (($config{'SEE'}) && ($config{'SEE'} =~ /false/i)) { $config{'SEE'} = 0; }

	$config{'operation_mode'} = $properties_obj->get_value("Parameters", "MODE") or die "Error, need MODE";
	die "Error, can't parse mode $config{'operation_mode'}\n" unless ($config{'operation_mode'} =~ /RELATIVE_POSITION|GENOME_COORDINATE/);

	$config{'inputFile'} = $properties_obj->get_value("MatchPairs", "Data");
	$config{'matchFormat'} = $properties_obj->get_value("MatchPairs", "Format");

	$config{'GAP_LENGTH'} = $properties_obj->get_value("Parameters", "GAP_LENGTH") or die "Error, need GAP_LENGTH";
	$config{'GAP_OPEN_PENALTY'} = $properties_obj->get_value("Parameters", "GAP_OPEN");
	die "Error, need GAP_OPEN" unless ((defined $config{'GAP_OPEN_PENALTY'}) && ($config{'GAP_OPEN_PENALTY'} <= 0));
	$config{'GAP_EXTENSION_PENALTY'} = $properties_obj->get_value("Parameters", "GAP_EXTEND");
	die "Error, need GAP_EXTEND" unless ((defined $config{'GAP_EXTENSION_PENALTY'}) && ($config{'GAP_EXTENSION_PENALTY'} <= 0));

	$config{'MAX_MATCH_SCORE'} = $properties_obj->get_value("Parameters", "MAX_MATCH_SCORE");
	$config{'MAX_DIST_BETWEEN_MATCHES'} = $properties_obj->get_value("Parameters", "MAX_DIST_BETWEEN_SYN_PAIRS") or die "Error, need MAX_DIST_BETWEEN_SYN_PAIRS";

	$config{'INCLUDE_SELF'} = $properties_obj->get_value("Parameters", "INCLUDE_SELF_COMPARISONS") || 0;
	if ($config{'INCLUDE_SELF'} =~ /false/) { $config{'INCLUDE_SELF'} = 0; }

	$config{'CONSTANT_MATCH_SCORE'} = $properties_obj->get_value("Parameters", "CONSTANT_MATCH_SCORE");

	# for auto-determine min alignment score cutoff
	if ($config{'CONSTANT_MATCH_SCORE'}) { $config{'MAX_MATCH_SCORE'} = $config{'CONSTANT_MATCH_SCORE'}; }

	$config{'MIN_NUM_ALIGNED_PAIRS'} = $properties_obj->get_value("Parameters", "MIN_ALIGNED_PAIRS") or die "Error, need MIN_ALIGNED_PAIRS";

	$config{'max_e_value'} = $properties_obj->get_value("Parameters", "MAX_EVALUE");
	die "Error, need MAX_EVALUE" unless (defined $config{'max_e_value'});

	$config{'MIN_ALIGNMENT_SCORE'} = $properties_obj->get_value("Parameters", "MIN_ALIGNMENT_SCORE");
	unless ($config{'MIN_ALIGNMENT_SCORE'}) {
		$config{'MIN_ALIGNMENT_SCORE'} = int ($config{'MIN_NUM_ALIGNED_PAIRS'} * 0.5 * $config{'MAX_MATCH_SCORE'});
	}

	$config{'NOISE_FILTER_DIST'} = $properties_obj->get_value("NoiseFilter", "BEST_MATCH_AGGREGATE_DIST");

	my %organism_to_gff_file;
	my @organisms = $properties_obj->get_section_attributes("GeneAnnotations");
	foreach my $organism (@organisms) {
		my $gff_file = $properties_obj->get_value("GeneAnnotations", $organism) or die "Error, cannot identify gff file for $organism";
		$organism_to_gff_file{$organism} = $gff_file;
	}
	$config{'orthologs_file'} = $properties_obj->get_value("Orthologs", "Data");
	return (\%config, \%organism_to_gff_file);
}
