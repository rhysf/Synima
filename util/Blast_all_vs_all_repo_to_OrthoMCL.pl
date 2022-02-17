#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use RBH_manager;
use File::Basename;
use Cwd;
use Synima;
use read_GFF;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <Repo_spec>\n
Optional: -t Type (PEP/CDS) [PEP]
          -o Out directory [OMCL_outdir]\n";
our($opt_r, $opt_t, $opt_o);
getopt('rtofsdm');
die $usage unless($opt_r);
if(!defined $opt_t) { $opt_t = 'PEP'; }
if(!defined $opt_o) { $opt_o = 'OMCL_outdir'; }
die "-t needs to be PEP or CDS: $opt_t\n" unless ($opt_t =~ m/^(PEP|CDS)$/);

# Scripts needed
my $m8_to_orthomcl = "$Bin/support_scripts/Blast_m8_to_OrthoMCL_gg_and_bpo_input.pl";
my $OrthoMCL = "$Bin/support_scripts/OrthoMCL.pl";
foreach($m8_to_orthomcl, $OrthoMCL){ die "Cannot find $_ : $!\n" unless(-e $_); }

# Input
my $all_annotations = "$opt_r.all.GFF3";

# Output directory and files
`mkdir $opt_o`;
my $all_blast_pairs = "$opt_o/all_blast_pairs.m8";
my $genome_codes_out = "$opt_o/for_omcl.genome_codes";
my $blast_pairs_gcoded = "$opt_o/for_omcl.Gcoded.m8";

# Get blast data into a single file
&retrieve_blast_pairs($opt_r, $opt_t, $all_blast_pairs);

# Assign genome codes to genes for omcl
my ($trans_id_to_genome, $genome_to_code) = gfffile::save_genome_codes_from_synima_parsed_gff3($all_annotations);
$genome_to_code = &assign_genome_codes($genome_to_code, $genome_codes_out);
&write_gcoded_m8_and_sort($trans_id_to_genome, $genome_to_code, $all_blast_pairs, $blast_pairs_gcoded);

# Convert m8 format to orthomcl input format
my $cmd = "perl $m8_to_orthomcl $blast_pairs_gcoded 4 $opt_o/omcl_in";
synima::process_cmd($cmd);

# Run orthomcl
$cmd = "perl $OrthoMCL --mode 4 --bpo_file $opt_o/omcl_in.bpo --gg_file $opt_o/omcl_in.gg 2>&1 | tee $opt_o/omcl.log";
synima::process_cmd($cmd);

# Move results from OrthoMCL directory to CWD
my $omcl_log_text = `cat $opt_o/omcl.log`;
$omcl_log_text =~ /Final ORTHOMCL Result: (\S+)/ or die "Error, cannot extract final result filename from $opt_o/omcl.log";
my $omcl_results_file = $1 or die "Error, cannot extract final result filename from $opt_o/omcl.log";
my $base_filename = basename($omcl_results_file);
my $tmp_dirname = dirname($omcl_results_file);
my $new_filename = cwd() . "/$base_filename";
$cmd = "mv $omcl_results_file $new_filename";
synima::process_cmd($cmd);
warn "Relocated $omcl_results_file to $base_filename\n";
$cmd = "rm -r $tmp_dirname";
synima::process_cmd($cmd);
print "Done.\n";

sub retrieve_blast_pairs {
	my ($repo_spec, $type, $output) = @_; 

	# Genomes in lexical order from repo spec
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();
	@genomes = sort @genomes; 

	# Make each protein file blastable
	warn "retrieve_blast_pairs: Printing to $output\n";
	my @blast_cmds;
	for (my $i = 0; $i <= $#genomes; $i++) {
		my $genomeA = $genomes[$i];
		GENOMES2: for (my $j = $i; $j <= $#genomes; $j++) {
			my $genomeB = $genomes[$j];
			my $genomeA_vs_genomeB_blast_file = &RBH_manager::get_blast_output_file($data_manager, $genomeA, $genomeB, $type);
			die "Error, no required blast output file: $genomeA_vs_genomeB_blast_file\n" if (! -s $genomeA_vs_genomeB_blast_file);

			# Print blast commands (only for OMCL, and not RBH)
			&retrieve_hits($genomeA_vs_genomeB_blast_file, $output);

			# Just do it once (OMCL)
			next GENOMES2 if ($genomeA eq $genomeB);
 
			my $genomeB_vs_genomeA_blast_file = &RBH_manager::get_blast_output_file($data_manager, $genomeB, $genomeA, $type);
			die "Error, no required blast output file: $genomeB_vs_genomeA_blast_file\n" if (! -s $genomeB_vs_genomeA_blast_file);
			&retrieve_hits($genomeB_vs_genomeA_blast_file, $output);
		}
	}
	return;
}

sub retrieve_hits {
	my ($file, $output) = @_;
	open my $fh, '<', $file or die "Error, cannot open file $file : $!\n";
	open my $ofh, '>>', $output or die "Error, cannot open file $output : $!\n";
	while(my $line=<$fh>) {
		# Error check
		my @bits = split /\t/, $line;
		die "$file is incomplete/corrupted for line consisting of: $line\nRemake file or delete line.\n" if(scalar(@bits) ne 12);
		# Print
		print $ofh $line;
	}
	close $fh;
	close $ofh;
	return;
}

sub assign_genome_codes {
	my ($genome_to_code, $outfile) = @_;
	my $counter = 0;
	warn "Assigning genome codes...\n";
	open my $ofh, '>', $outfile or die "Error, cannot write to $outfile : $!\n";
	foreach my $genome (keys %{$genome_to_code}) {
		$counter++;
		my $code = sprintf("G%03x", $counter);
		$$genome_to_code{$genome} = $code;
		warn "$genome => $code\n";
		print $ofh "$genome\t$code\n";
	}
	close $ofh;
	die "assign_genome_codes failed. Rerun.\n" if(scalar(keys(%{$genome_to_code})) < 1);
	return $genome_to_code;
}

sub write_gcoded_m8_and_sort {
	my ($trans_id_to_genome, $genome_to_code, $blast_hits, $outfile) = @_;

	# Write Gcoded version
	open my $fh, '<', $blast_hits or die "Cannot open $blast_hits : $!\n";
	open my $ofh, '>', "$outfile.tmp" or die "Error, cannot write to file $outfile.tmp : $!\n";
	warn "Writing Gcoded version of blast.m8...\n";
	while (my $line=<$fh>) {
		chomp $line;
		my @x = split /\t/, $line;
		my ($accA, $accB) = @x;
		my $genomeA = $$trans_id_to_genome{$accA} or die "Error, no genome for $accA saved from all_annotations.gff3 (check settings and rerun): $!\n";
		my $genomeB = $$trans_id_to_genome{$accB} or die "Error, no genome for $accB saved from all_annotations.gff3 (check settings and rerun): $!\n";
		my $codeA = $$genome_to_code{$genomeA};
		my $codeB = $$genome_to_code{$genomeB};
		$x[0] = "$codeA|$accA";
		$x[1] = "$codeB|$accB";
		print $ofh join("\t", @x) . "\n";
	}
	close $fh;
	close $ofh;

	# Sort the Gcoded version
	my $cmd = "sort -T . -S 2G -k1,1 -k12,12gr $outfile.tmp > $outfile";
	synima::process_cmd($cmd);
	unlink("$outfile.tmp");
	warn "Done\nSee file: $outfile\n";
}
