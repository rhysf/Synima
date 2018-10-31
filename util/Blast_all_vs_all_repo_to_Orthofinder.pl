#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use DataSpecFileParser;
use RBH_manager;
use File::Basename;
use Data::Dumper;
use Cwd;
use Synima;
use read_GFF;
use read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <Repo_spec>\n
Optional: -t Type (PEP/CDS) [PEP]
          -o Out directory [Orthofinder_outdir]\n";
our($opt_r, $opt_t, $opt_o);
getopt('rto');
die $usage unless($opt_r);
if(!defined $opt_t) { $opt_t = 'PEP'; }
if(!defined $opt_o) { $opt_o = 'Orthofinder_outdir'; }
die "-t needs to be PEP or CDS: $opt_t\n" unless ($opt_t =~ m/^(PEP|CDS)$/);

# Dependencies
my $uname = `uname`;
chomp $uname;
my $fastme = "$Bin/support_scripts/fastme.$uname";
my $orthofinder = "$Bin/support_scripts/orthofinder.$uname";
my $mcl = "$Bin/support_scripts/mcl.$uname";
#my $orthofinder = "$Bin/support_scripts/orthofinder.py";
foreach($fastme, $mcl, $orthofinder) { die "Cannot find $_ : $!\n" unless(-e $_); }
if(-e "$Bin/support_scripts/mcl") { `rm $Bin/support_scripts/mcl`; }
if(-e "$Bin/support_scripts/fastme") { `rm $Bin/support_scripts/fastme`; }
`ln -s $mcl $Bin/support_scripts/mcl`;
`ln -s $fastme $Bin/support_scripts/fastme`;
local $ENV{PATH} = "$ENV{PATH}:$Bin/support_scripts/";

# Parse repo spec into struct
warn "Parsing repo spec...\n";
my $data_manager = new DataSpecFileParser($opt_r);

# Move and convert files
die "Error: Directory $opt_o already exists. Delete and relaunch\n" if(-d $opt_o);
`mkdir $opt_o`;

# Assign genome codes to genes for orthofinder
my $genome_codes_out = "$opt_o/SpeciesIDs.txt";
my $genome_to_code = &assign_genome_codes_Orthofinder($data_manager, $opt_t, $genome_codes_out);

# Assign sequence IDs and save FASTA's
my $seq_IDs_to_new_IDs = &assign_sequence_IDs_and_save_FASTA($data_manager, $opt_t, $genome_to_code, $opt_o);

# Reformat BLAST reports
my $reformat_BLAST_reports_for_orthofinder = &retrieve_blast_pairs_orthofinder($opt_r, $opt_t, $genome_to_code, $seq_IDs_to_new_IDs, $opt_o);

# Run Orthofinder
warn "Running Orthofinder...\n";
my $run_of_cmd1 = `$orthofinder -b $opt_o`;
warn "$run_of_cmd1\n";
warn "Finished with: $?\n";

sub assign_sequence_IDs_and_save_FASTA {
	my ($data_manager, $dna_or_protein, $genome_to_code, $outfolder) = @_;
	warn "Assigning sequence IDs and save FASTA...\n";
	my @genomes = $data_manager->get_genome_list();
	my $sequence_ids = "$outfolder/SequenceIDs.txt";
	my %seq_ids_to_new_ids;
	open my $ofh, '>', $sequence_ids or die "Error, cannot write to $sequence_ids : $!\n";
	for(my $i=0; $i<scalar(@genomes); $i++) {
		my $genome = $genomes[$i];
		my $fasta = $data_manager->get_data_dump_filename($genome, $dna_or_protein);
		$fasta .= ('.synima-parsed.' . $dna_or_protein);
		die "Cannot find $fasta : $!\n" if(! -e $fasta);
		my $code = $$genome_to_code{$genome};
		die "Cannot find a code for $genome : $!\n" if(!defined $$genome_to_code{$genome});
		#warn "do stuff with $fasta and $code?\n";

		# Save FASTA
		my $seq_counter = 0;
		my ($sequences, $descriptions, $order) = fastafile::fasta_id_to_seq_hash($fasta);
		my $outfile2 = "$outfolder/Species$code.fa";
		open my $ofh2, '>', $outfile2 or die "Cannot open outfile $outfile2 : $!\n";
		foreach my $id(keys %{$sequences}) {
			my $seq = $$sequences{$id};
			my $new_id = "$code\_$seq_counter";
			#warn "new id = $new_id and seq = $seq\n";
			print $ofh "$new_id: $id\n";
			die "already found a gene with id $id in another genome. Make unique IDs for each genome.\n" if(defined $seq_ids_to_new_ids{$id});
			$seq_ids_to_new_ids{$id} = $new_id;
			print $ofh2 ">$new_id\n$seq\n";
			$seq_counter++;
		}
		close $ofh2;
	}
	close $ofh;
	return \%seq_ids_to_new_ids;
}

sub assign_genome_codes_Orthofinder {
	my ($data_manager, $dna_or_protein, $outfile) = @_;

	my %genome_to_code;
	my $counter = 0;
	warn "Assigning genome codes...\n";
	my @genomes = $data_manager->get_genome_list();
	open my $ofh, '>', $outfile or die "Error, cannot write to $outfile : $!\n";
	for(my $i=0; $i<scalar(@genomes); $i++) {
		my $genome = $genomes[$i];
		my $fasta = $data_manager->get_data_dump_filename($genome, $dna_or_protein);
		$fasta .= ('.synima-parsed.' . $dna_or_protein);
		die "Cannot find $fasta : $!\n" if(! -e $fasta);
		my $base_fasta = basename($fasta);
		$genome_to_code{$genome} = $counter;
		#warn "$counter => $fasta\n";
		print $ofh "$counter: $base_fasta\n";
		$counter++;
	}
	close $ofh;
	die "assign_genome_codes failed. Rerun.\n" if(scalar(keys(%genome_to_code)) < 1);
	return \%genome_to_code;
}

sub retrieve_blast_pairs_orthofinder {
	my ($repo_spec, $type, $genome_to_code, $seq_ids_to_new_ids, $outfolder) = @_; 
	warn "retrieve_blast_pairs_orthofinder...\n";

	# Genomes in lexical order from repo spec
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();
	@genomes = sort @genomes; 

	# Make each protein file blastable
	my @blast_cmds;
	for (my $i = 0; $i <= $#genomes; $i++) {
		my $genomeA = $genomes[$i];
		my $genomeA_code = $$genome_to_code{$genomeA};
		GENOMES2: for (my $j = $i; $j <= $#genomes; $j++) {
			my $genomeB = $genomes[$j];
			my $genomeB_code = $$genome_to_code{$genomeB};
			my $genomeA_vs_genomeB_blast_file = &RBH_manager::get_blast_output_file($data_manager, $genomeA, $genomeB, $type);
			die "Error, no required blast output file: $genomeA_vs_genomeB_blast_file\n" if (! -s $genomeA_vs_genomeB_blast_file);

			# Print blast commands (only for orthofinder)
			#die "genomeA $genomeA $genomeA_code -> genomeB $genomeB $genomeB_code\n";
			my $output = ($outfolder . "/Blast" . $genomeA_code . "_" . $genomeB_code . ".txt");
			&retrieve_hits_orthofinder($genomeA_vs_genomeB_blast_file, $output, $seq_ids_to_new_ids);

			# Just do it once (Orthofinder)
			next GENOMES2 if ($genomeA eq $genomeB);
 
			my $genomeB_vs_genomeA_blast_file = &RBH_manager::get_blast_output_file($data_manager, $genomeB, $genomeA, $type);
			die "Error, no required blast output file: $genomeB_vs_genomeA_blast_file\n" if (! -s $genomeB_vs_genomeA_blast_file);
			$output = ($outfolder . "/Blast" . $genomeB_code . "_" . $genomeA_code . ".txt");
			&retrieve_hits_orthofinder($genomeB_vs_genomeA_blast_file, $output, $seq_ids_to_new_ids);
		}
	}
	return;
}

sub retrieve_hits_orthofinder {
	my ($file, $output, $seq_ids_to_new_ids) = @_;
	warn "retrieve_hits_orthofinder: $file -> $output\n";
	open my $fh, '<', $file or die "Error, cannot open file $file : $!\n";
	open my $ofh, '>>', $output or die "Error, cannot open file $output : $!\n";
	while(my $line=<$fh>) {
		# Error check
		my @bits = split /\t/, $line;
		die "$file is incomplete/corrupted for line consisting of: $line\nRemake file or delete line.\n" if(scalar(@bits) ne 12);

		# Change IDs
		my $new_id1 = $$seq_ids_to_new_ids{$bits[0]};
		my $new_id2 = $$seq_ids_to_new_ids{$bits[1]};
		die "Cannot find assigned IDs from SequenceIDs.txt for $line ($new_id1 and $new_id2)\n" if((!defined $new_id1) || (!defined $new_id2));
		$bits[0] = $new_id1;
		$bits[1] = $new_id2;
		my $new_line = join "\t", @bits;
		# Print
		print $ofh $new_line;
		#print $new_line;
		#die "end here! check output file\n";
	}
	close $fh;
	close $ofh;
	return;
}
