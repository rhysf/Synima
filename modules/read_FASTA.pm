package fastafile;
use strict;
use Bio::SeqIO;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);

### rfarrer@broadinstitute.org

sub fasta_to_struct {
	my $input = $_[0];
	my %struct;
	$struct{'filename'} = $input;
	warn "fasta_to_struct: saving from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		my $seq = $seq_obj->seq;
		my $desc = $seq_obj->description;

		# Save
		$struct{'seq'}{$id} = $seq;
		$struct{'desc'}{$id} = $desc;
		push @{$struct{'order'}}, $id;
	}
	return \%struct;
}

sub fasta_id_to_seq_length_hash {
	my $input = $_[0];
	my (%lengths);
	warn "fasta_id_to_seq_length_hash: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
    	my $seq = $seq_obj->seq;
		my $length = length($seq);
		$lengths{$id} = $length;
	}
	return (\%lengths);
}

sub fasta_to_total_seq_length {
	my $input = $_[0];
	my ($lengths) = 0;
	warn "fasta_to_total_seq_length: saving sequences from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
    	my $seq = $seq_obj->seq;
		$lengths += length($seq);
	}
	return ($lengths);
}

sub fasta_id_to_order_array {
	my $input = $_[0];
	my @order;
	warn "fasta_id_to_order_array: saving order from $input...\n";
	my $inseq = Bio::SeqIO->new('-file' => "<$input",'-format' => 'fasta');
	while (my $seq_obj = $inseq->next_seq) { 
		my $id = $seq_obj->id;
		push @order, $id;
	}
	return (\@order);
}

sub fasta_broad_format_to_struct {
	my $file = $_[0];

	my %all_proteins;
	warn "fasta_broad_format_to_struct: $file...\n";
	open my $fh, '<', $file or die "Unable to open file $file : $!\n";
	PROTEIN: while(my $line=<$fh>){
		chomp $line;
		next PROTEIN if($line !~ m/>/);
		my ($transcript_id, $gene_id, $locus_name, $func_annot, $genome, $analysis) = &parse_protein_file_line($line);
		die "genome not defined : $genome\n" unless($genome);
		$all_proteins{$transcript_id} = { transcript_id => $transcript_id, gene_id => $gene_id, locus_name => $locus_name, func_annot => $func_annot, genome => $genome, analysis => $analysis };
	}
   	close $fh; 
	return \%all_proteins;
}

sub create_transcript_functinal_annotation_map {
	my ($file_in, $file_out) = @_;
	open my $fh, '<', $file_in or die "Unable to open file $file_in : $!\n";
	open my $ofh, '>', $file_out or die "Unable to open file $file_out : $!\n";
	warn "create_transcript_functinal_annotation_map: $file_out\n";
	while(my $line=<$fh>) {
		chomp $line;
		next if ($line !~ m/^>/);
		my ($transcript_id, $gene_id, $locus_name, $func_annot, $genome, $analysis) = &parse_protein_file_line($line);
		print $ofh "$transcript_id\t$gene_id\t$locus_name\t$func_annot\n";
	}
	close $fh;
	close $ofh;
	return 1;
}

sub split_fasta_seq_dictionary_by_species {
	my ($repo_spec, $type) = @_;

	# Parse repo spec file
	my $data_manager = new DataSpecFileParser($repo_spec);
	my @genomes = $data_manager->get_genome_list();

	# Input seq dictionary
	my $input = "$repo_spec.all.$type";
	die "Cannot find $input. Re-run Create_full_repo_sequence_databases.pl\n" if(! -e $input);
	my $fasta = fastafile::fasta_to_struct($input);

	# Split seq dictionary by genomes
	warn "split_fasta_seq_dictionary_by_species: split $input\n";
	foreach my $genome(@genomes) {

		# Outfile
		my $count = 0;
		my $output_file = $data_manager->get_data_dump_filename($genome, $type);
		$output_file .= ".synima-parsed.$type";
		warn "split_fasta_seq_dictionary_by_species: splitting $genome -> $output_file...\n";

		# Print synima-parsed.CDS/PEP files with updated ids
		open my $ofh, '>', $output_file or die "Cannot open output file $output_file : $!\n";
		foreach my $transcript_id(keys %{$$fasta{'seq'}}) {

			# gene/protein ID
			my $id_line = ">$transcript_id";
			if(defined $$fasta{'desc'}{$transcript_id}) { $id_line .= " $$fasta{'desc'}{$transcript_id}"; }

			# gene/protein sequence
			my ($transcript_id, $gene_id, $locus_name, $func_annot, $genome_found, $analysis) = &parse_protein_file_line($id_line);
			next if($genome_found ne $genome);
			$count++;
			print $ofh ">$transcript_id\n$$fasta{'seq'}{$transcript_id}\n";
		}
		close $ofh;
		warn "split_fasta_seq_dictionary_by_species: $count printed for $genome\n";
	}
	return 1;
}

### local sub routines

#Broad Format = >7000011728610201 gene_id=7000011728610200 locus=None name="flagellum-specific ATP synthase" genome=Esch_coli_MGH121_V1 analysisRun=Esch_coli_MGH121_V1_POSTPRODIGAL_2	
sub parse_protein_file_line {
	my $line =  shift @_;
	my @cols = split " ", $line;

	my $transcript_id = shift @cols;	 	
	die "Error:Unable to parse line $line from proteins file. Transcript_id not defined.\n" if(!defined $transcript_id);
	$transcript_id =~ s/^>//;

	my $gene_id = shift @cols;
	die "Error:Unable to parse line $line from proteins file. Gene_id not defined.\n" if(!defined $gene_id);
	$gene_id =~ s/gene_id=//;
	
	my $locus_name = shift @cols;
	die "Error:Unable to parse line $line from proteins file. Locus name not defined.\n" if(!defined $locus_name);
	$locus_name =~ s/locus=//;

	my $func_annot = '';
	my $index = 0;
	while($index < scalar(@cols)) {
		$func_annot .= $cols[$index] . " ";
		$index++;
		#warn "index = $index, func_annot = $func_annot\n";
		#last if($func_annot =~ /^name=""/);
		last if($func_annot =~ /^name="[\w\W]+" $/);
	}

	# If we got to the end of the line...
	die "parse_protein_file_line: Error: Unable to parse line $line from proteins file\n" if ($index - 1 ) == scalar(@cols);

	$func_annot =~ s/name="//g;
	$func_annot =~ s/" //g;

	my $genome = $cols[$index];
	die "Unable to find genome in $line at index $index\n" if(!defined $genome);
	$index++;	
	$genome =~ s/genome=//;

	my $analysis = $cols[$index];
	$analysis =~ s/analysisRun=//;
	return ($transcript_id, $gene_id, $locus_name, $func_annot, $genome, $analysis);
}

1;
