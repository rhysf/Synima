#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use read_FASTA;
use Run_Bsub;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling);

### r.farrer@exeter.ac.uk

# Opening statements
my $usage = "Usage: perl $0 -r <repo> -c <gene_clusters orthologAnalysisRun.clusters (ortholog data dump output)>\n
Optional: -d use codon alignments instead of protein alignments [0]
          -g all computes are executed on the computing grid (y/n) [n]
          -s file containing the list of genomes to restrict alignments and trees to []
          -t Make trees for every orthogroup (y/n) [n]\n";
our($opt_r, $opt_c, $opt_d, $opt_g, $opt_s, $opt_t);
getopt('rcdgst');
die $usage unless($opt_r && $opt_c);
if(!defined $opt_d) { $opt_d = 0; }
if(!defined $opt_g) { $opt_g = 'n'; }
if(!defined $opt_t) { $opt_t = 'n'; }

# Programs
my $run_cmds_on_grid = "$FindBin::Bin/support_scripts/Run_cmds_on_grid.py";
my $peptide_to_malign  = "$FindBin::Bin/support_scripts/FASTA_peptide_to_malign.pl";
foreach($peptide_to_malign) { die "Cannot find $_ : $!\n" unless (-e $_); } 
if($opt_g ne 'n') { 
	foreach($run_cmds_on_grid) { die "Cannot find $_ : $!\n" unless (-e $_); } 
}

# Check that MUSCLE is installed
my $muscle_path = `which muscle 2>/dev/null`;
chomp($muscle_path);

if (!$muscle_path || !-x $muscle_path) {
    die "Error: MUSCLE not found in PATH or not executable.\n";
}

# Check version (muscle -version is available only in v5)
my $version_output = `muscle -version 2>&1`;
if ($version_output =~ /muscle\s+v?([\d\.]+)/i) {
    my $version = $1;
    $version =~ s/[^0-9.]//g;
    $version =~ s/\.$//;
    print "Found MUSCLE version: $version\n";

    if ($version < 5) {
        die "Error: MUSCLE v5+ is required (you have v$version).\n" .
            "Please install MUSCLE v5: https://github.com/rcedgar/muscle\n";
    }
} else {
    die "Error: Could not determine MUSCLE version. MUSCLE v5 is required.\n" .
        "Try running 'muscle -version' manually.\n";
}

# Out directory
my $malign_dir = basename($opt_c) . ".MALIGN_DIR";
if ($opt_s) { $malign_dir = basename($opt_s) . ".MALIGN_DIR"; }

my %proteins_by_genome; # all protein sequences
my %CDS_by_genome; # only if doing codon alignments
my $data_manager = new DataSpecFileParser($opt_r);

# Perform data retrievals
my @genomes;
if ($opt_s) {
	@genomes = `cat $opt_s`;
	foreach my $genome (@genomes) { $genome =~ s/\s+//g; }
	@genomes = grep { /\w/ } @genomes; # ignore empty entries
}
else { @genomes = $data_manager->get_genome_list(); }

main: {

	# Go through each genome saving sequences
	foreach my $genome (@genomes) {

		# Save all the proteins
		my $protein_file_base = $data_manager->get_data_dump_filename($genome, 'PEP');
		my $protein_file = "$protein_file_base.synima-parsed.PEP";
		my $protein_seqs = fastafile::fasta_to_struct($protein_file);
		$proteins_by_genome{$genome} = $protein_seqs;

		# Save CDS as well
		if ($opt_d) {
			my $cds_file_base = $data_manager->get_data_dump_filename($genome, 'CDS');
			my $cds_file = "$cds_file_base.synima-parsed.PEP";
			my $cds_seqs = fastafile::fasta_to_struct($cds_file);
			$CDS_by_genome{$genome} = $cds_seqs;
		}
	}
	my %clusters = &parse_gene_clusters($opt_c, \@genomes);


	# Make the MALIGN output directory
	if (! -d $malign_dir) {
		mkdir($malign_dir) or die "Error, cannot mkdir $malign_dir : $!\n";
	}

	# For each ortholog cluster, make a protein and CDS FASTA.
	my @pep_files;
	foreach my $cluster_id (keys %clusters) {
		warn "// processing Cluster: $cluster_id\n";

		# Write peptide multiple alignment fasta file
		my $outfile = "$malign_dir/$cluster_id.pep";
		if (-s $outfile) { warn "-$outfile already exists, not overwriting it\n"; }
		else { 
			open my $ofh, '>', $outfile or die "Error, cannot write to file $outfile : $!\n";
			foreach my $genome (sort keys %{$clusters{$cluster_id}}) {
				my @trans = @{$clusters{$cluster_id}->{$genome}};
				foreach my $trans_id (@trans) {
		  			#my $pep = uc $proteins_by_genome{$genome}->{$trans_id} or die "Error, no protein sequence for $genome, $trans_id\n";
		  			my $pep = uc $proteins_by_genome{$genome}{'seq'}{$trans_id} or die "Error, no protein sequence for $genome, $trans_id\n";
		  			print $ofh ">$trans_id\n$pep\n";
		   		}
	   		}
	   		close $ofh;
		}

		if ($opt_d) {
	   		my $cds_outfile = "$malign_dir/$cluster_id.cds";
	   		if (-s $cds_outfile) { warn "-$cds_outfile already exists, not overwriting it"; } 
			else {
				my $cds_text = "";
				foreach my $genome (sort keys %{$clusters{$cluster_id}}) {
					my @trans = @{$clusters{$cluster_id}->{$genome}};
					foreach my $trans_id (@trans) {
			  
						my $cds = uc $CDS_by_genome{$genome}{'seq'}{$trans_id} or die "Error, no cds sequence for $genome, $trans_id\n";
						$cds_text .= ">$trans_id\n$cds\n";
					}
				}
				open my $ofh, '>', $cds_outfile or die "Error, cannot write to $cds_outfile\n";
				print $ofh $cds_text;
				close $ofh;
			}
		}
		push (@pep_files, $outfile);
	}

	# Generate the multiple alignments
	my @seqs_to_align = @pep_files;
	my @cluster_cmds;
	foreach my $seqs_to_align_file (@seqs_to_align) {

		my $cmd = "perl $peptide_to_malign -s $seqs_to_align_file -t $opt_t";

		# make codon alignment instead
		if ($opt_d) { $cmd .= " -c 1"; }

		# run on grid
		if ($opt_g ne 'n') {
			push (@cluster_cmds, $cmd);
		}
		else { &process_cmd($cmd); }
	}

	# run on grid
	if ($opt_g ne 'n') {
		my $cmd_outfile = "$malign_dir/cluster_cmds";
		open (my $ofh, ">$cmd_outfile") or die "Error, cannot write to $cmd_outfile";
		foreach my $ind_cmd (@cluster_cmds) {
			print $ofh "$ind_cmd\n";
		}
		close $ofh;

		my $dotkits_file = "$malign_dir/dotkits_list";
		open (my $dfh, ">$dotkits_file") or die "Error, cannot write to $dotkits_file";
		print $dfh "UGER\nPython-2.7\n";
		close $dfh;

		warn "-running protein multiple alignment and tree building commands\n";
		my $run_grid_cmd = "$run_cmds_on_grid --platform UGER --queue long --mem 4 --throttle_nodes 99 --cmds_per_node 1 --dotkits $dotkits_file $cmd_outfile";
		&process_cmd($run_grid_cmd);
	}
	print "Done.\n\n";
	exit(0);
}

sub process_cmd {
	my ($cmd) = @_;
	warn "CMD: $cmd\n";
	my $ret = system($cmd);
	if ($ret) {
		warn "Error, cmd: $cmd died with ret: $ret";	
	}
	return;
}

sub parse_gene_clusters {
	my ($gene_clusters_file, $genomes_aref) = @_;
	
	my %genomes = map { + $_ => 1 } @$genomes_aref;
	my %cluster_to_data;
	
	open (my $fh, $gene_clusters_file) or die "Error, cannot open file $gene_clusters_file";
	CLUSTERS: while (my $line=<$fh>) {
		next CLUSTERS unless ($line =~ m/\w/);
		my @x = split /\t/, $line;
		my ($cluster_id, $genome, $gene_set_name, $trans_id, $gene_id, $locus, @rest) = @x;

		# ignore those entries that are not in the list
		next CLUSTERS unless ($genomes{$genome});
		push (@{$cluster_to_data{$cluster_id}->{$genome}}, $trans_id);
	}
	return(%cluster_to_data);
}
