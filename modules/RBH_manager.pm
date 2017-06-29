package RBH_manager;
use strict;
use warnings;
use Carp;
use DataSpecFileParser;

sub get_blast_output_directory {
	my ($data_manager, $genomeA, $genomeB, $type) = @_;
	my $genomeA_repo_dir = $data_manager->get_genome_repo_dir($genomeA);
	my $blast_out_dir = "$genomeA_repo_dir/RBH_blast_$type";
	return($blast_out_dir);
}

sub get_blast_output_file {
	my ($data_manager, $genomeA, $genomeB, $type) = @_;
	my $blast_output_dir = &get_blast_output_directory($data_manager, $genomeA, $genomeB, $type);
	my $outfile = "$blast_output_dir/vs_$genomeB.blast.m8";
	return($outfile);
}

sub get_RBH_output_file {
	my ($data_manager, $genomeA, $genomeB, $type) = @_;
	($genomeA, $genomeB) = sort ($genomeA, $genomeB);
	my $blast_output_file = &get_blast_output_file($data_manager, $genomeA, $genomeB, $type);
	my $rbh_file = "$blast_output_file.RBH";
	return($rbh_file);
}

1; #EOM
