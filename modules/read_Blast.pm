package blastfile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use File::Which;
use FindBin qw($Bin);
use lib "$Bin";
use DataSpecFileParser;
use RBH_manager;
use Synima;

### rfarrer@broadinstitute.org

sub make_all_vs_all_blast_search_cmds_from_repo {
	my ($repo_spec, $type, $num_matches, $num_hits, $evalue, $outfile) = @_;
	my $blast_prog = ($type eq "PEP") ? "blastp" : "blastn";

	# BLAST version
	my $BLAST_version = "legacy BLAST";
	my $formatdb_path = which 'formatdb';
	my $makeblastdb_path = which 'makeblastdb';
	if($formatdb_path eq '') { 
		if($makeblastdb_path eq '') { die "Legacy BLAST (formatdb) and BLAST+ (makeblastdb) not found. Make sure one is in your PATH\n"; } 
		else { $BLAST_version = "BLAST+"; }
	}
	warn "Using $BLAST_version...\n";

	# Parse repo spec file
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();

	# Make commands for all-vs-all-search
	open my $blast_cmds_ofh, '>', $outfile or die "Error, cannot write to $outfile\n";
	foreach my $genomeA (@genomes) {
		my $genomeA_file = $data_manager->get_data_dump_filename($genomeA, $type);
		$genomeA_file .= ".synima-parsed.$type";
		foreach my $genomeB (@genomes) {
			my $genomeB_file = $data_manager->get_data_dump_filename($genomeB, $type);
			$genomeB_file .= ".synima-parsed.$type";

			# Make blastable
			&make_genome_blastable($genomeB_file, $type, $BLAST_version);

			# Make output directory and specify outfile
			my ($genomeA_vs_genomeB_blast_file, $ignore) = &make_outdir_and_output_for_blast($data_manager, $genomeA, $genomeB, $type);
			next if($ignore eq 1);

			# Make Blast command
			my $num_hits = ($genomeA eq $genomeB) ? $num_hits : $num_matches;
			my $cmd;
			if($BLAST_version eq 'legacy BLAST') { $cmd = "blastall -p $blast_prog -i $genomeA_file -d $genomeB_file -m 8 -v $num_hits -b $num_hits -e $evalue > $genomeA_vs_genomeB_blast_file"; }
			else { $cmd = "$blast_prog -query $genomeA_file -db $genomeB_file -outfmt 6 -max_target_seqs $num_hits -evalue $evalue > $genomeA_vs_genomeB_blast_file" }
			print $blast_cmds_ofh "$cmd\n";
		}
	}
	close $blast_cmds_ofh;
	return 1;
}

sub make_genome_blastable {
	my ($genome_file, $type, $BLAST_version) = @_;
	my $blastable_index;
	if($type eq 'PEP') { $blastable_index = "$genome_file.pin"; }
	else { $blastable_index = "$genome_file.nin"; }
	unless (-s $blastable_index) {
		if($BLAST_version eq 'legacy BLAST') {
			my $blast_cmd = "formatdb -i $genome_file -p ";
			if($type eq 'PEP') { $blast_cmd .= 'T'; }
			else { $blast_cmd .= 'F'; }
			synima::process_cmd($blast_cmd);
		} else {
			my $blast_cmd = "makeblastdb -in $genome_file -dbtype ";
			if($type eq 'PEP') { $blast_cmd .= 'prot'; }
			else { $blast_cmd .= 'nucl'; }
			synima::process_cmd($blast_cmd);
		}
	}
	return;
}

sub make_outdir_and_output_for_blast {
	my ($dm, $genomeA, $genomeB, $type) = @_;
	my $ignore = 0;

	# Make output directory
	my $genomeA_vs_genomeB_directory = RBH_manager::get_blast_output_directory($dm, $genomeA, $genomeB, $type);
	unless (-d $genomeA_vs_genomeB_directory) {
		mkdir ($genomeA_vs_genomeB_directory, 0775) or die "Error, cannot mkdir $genomeA_vs_genomeB_directory : $!\n";
	}

	# Do not over-write BLAST command if outfile already exists
	my $genomeA_vs_genomeB_blast_file = RBH_manager::get_blast_output_file($dm, $genomeA, $genomeB, $type);
	if (-s $genomeA_vs_genomeB_blast_file) {
		warn "-warning, $genomeA_vs_genomeB_blast_file already exists. Not over-writing it.\n";
		$ignore = 1;
	}

	return ($genomeA_vs_genomeB_blast_file, $ignore);
}

1;
