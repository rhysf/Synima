package synima;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use Cwd;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin";
use read_FASTA;
use read_DAG;
my $cwd = getcwd;

### rfarrer@broadinstitute.org

sub write_config {
	my @params = @_;
	my $aligncoords         = $params[0];  # a
	my $aligncoords_spans1  = $params[1];  # b
	my $config_file         = $params[2];  # c
	my $genome_extensions   = $params[3];  # e
	my $plot_contig_synteny = $params[4];  # g
	my $width_pixels        = $params[5];  # i
	my $height_pixels       = $params[6];  # j
	my $gene_list1          = $params[7];  # k
	my $gene_list2          = $params[8];  # l
	my $genome_labels       = $params[9];  # n
	my $gene_list3          = $params[10]; # o
	my $run_synima          = $params[11]; # r
	my $aligncoords_spans2  = $params[12]; # t
	my $aligncoords_spans3  = $params[13]; # u
	my $verbose             = $params[14]; # v
	my $left_margin         = $params[15]; # w
	my $order_of_genomes    = $params[16]; # x
	my $plot_ind_genes      = $params[17]; # z

	# Check input files can be read
	foreach($aligncoords, $aligncoords_spans1) { die "Cannot open file: $_ : $!\n" if(! -e $_); }

	# Make directory if needed
	my $dirname = dirname($config_file);
	if(! -d $dirname) { 
		warn "making output directory $dirname...\n";
		`mkdir $dirname`;
	}
	`touch $config_file`;

	# Find full paths for input files
	my $input_files = &check_input_file_make_config_entry($config_file, 'Self_location');
	$input_files .= &check_input_file_make_config_entry($aligncoords, 'aligncoords_location');
	$input_files .= &check_input_file_make_config_entry($aligncoords_spans1, 'aligncoords_spans_location1');
	$input_files .= &check_input_file_make_config_entry($aligncoords_spans2, 'aligncoords_spans_location2');
	$input_files .= &check_input_file_make_config_entry($aligncoords_spans3, 'aligncoords_spans_location3');
	$input_files .= &check_input_file_make_config_entry($gene_list1, 'gene_list1');
	$input_files .= &check_input_file_make_config_entry($gene_list2, 'gene_list2');
	$input_files .= &check_input_file_make_config_entry($gene_list3, 'gene_list3');
	$input_files =~ s/\n$//;

	# Save genome names, genome lengths, contig lengths, num of genomes, half way point and longest genome (xmax)
	my ($genome_names, $genome_lengths, $genome_contig_lengths, $num_genomes, $half_way_point_through_genomes, $xmax, $genome_contig_order) = &gather_genome_lengths_for_config($aligncoords_spans1, $genome_extensions, $verbose);

	# Genome labels
	if($genome_labels eq '') { $genome_labels = $genome_names; }
	my @genome_labels_split = split /,/, $genome_labels;
	my $number_of_genome_labels =  scalar(@genome_labels_split);
	die "Number of genomes ($num_genomes) not equal number of genome labels ($number_of_genome_labels): $genome_labels\n" if($number_of_genome_labels ne $num_genomes);

	# Height and width of figure
	if($height_pixels eq 'n') { 
		$height_pixels = ($num_genomes * 100);
		if($height_pixels < 400) { $height_pixels = 400; }
	}
	die "height in pixels not recognised as valid integer: $height_pixels\n" if($height_pixels !~ m/\d+/);
	die "width in pixels not recognised as valid integer: $width_pixels\n" if($width_pixels !~ m/\d+/);
	my $height_inches = sprintf("%.2f", ($height_pixels / 94));
	my $width_inches = sprintf("%.2f", ($width_pixels / 96));
	my $ymax = ($num_genomes - 1) + 0.2;

	# Make config file
	warn "Making config file $config_file...\n";
	open my $ofh, ">$config_file" or die "Cannot open $config_file : $!\n";
	print $ofh "### SynIma v1 config file
# Location of input and output files
$input_files
output_rscript=$config_file.Rscript
output_pdf=$config_file.pdf

# Genome info
genome_order=$genome_names
genome_labels=$genome_labels
genome_number_of=$num_genomes
genome_half_way_point=$half_way_point_through_genomes
$genome_contig_order

# Contig info
contig_min_length_for_label=400000

# Reverse compliment these contigs in the format reverse_compliment=genome1:chr1,genome1:chr2,genome2:chr1 etc.
reverse_compliment=

# Remove these letters from the begining of names in the format find1:replace1,find2:replace2
contig_name_replace=Supercontig:sc,Contig:c,Scaffold:,Chromosome:chr etc,

# Figure parameters
height_inches=$height_inches
width_inches=$width_inches
xmax=$xmax
ymax=$ymax
plot_synteny=$plot_contig_synteny
plot_individual_genes=$plot_ind_genes
left_margin=$left_margin
aligncoords_spans_colour1=azure4
aligncoords_spans_colour2=red
aligncoords_spans_colour3=blue
gene_list_colour1=cornflowerblue
gene_list_colour2=coral3
gene_list_colour3=darkcyan

### DONT EDIT BELOW THIS LINE
$genome_lengths
$genome_contig_lengths";
	close $ofh;

	# Now read what's written, and return
	my $synima_data = &read_config($config_file);
	return $synima_data;
}

sub read_config {
	my $config_file = $_[0];
	my %synima_data;
	open my $fh, '<', $config_file or die "Cannot open $config_file : $!\n";
	warn "Reading $config_file...\n";
	CONFIG: while(my $line=<$fh>) {
		chomp $line;
		next CONFIG if ($line =~ m/^#/);
		next CONFIG if ($line =~ m/^\n/);
		my @bits = split /=/, $line;
		next CONFIG if((!defined $bits[0]) || (!defined $bits[1]));

		# Save config
		if(scalar(@bits) eq 2) {
			#warn "$bits[0] -> $bits[1]\n";
			$synima_data{$bits[0]} = $bits[1]; 
		} else {
			#warn "$bits[0] -> $bits[1] -> $bits[2]\n";
			$synima_data{$bits[0]}{$bits[1]} = $bits[2]; 
		}
	}
	&config_sanity_check(\%synima_data);
	return \%synima_data;
}

sub config_sanity_check {
	my $synima_data = $_[0];
	foreach(qw(aligncoords_location aligncoords_spans_location1 output_rscript output_pdf genome_order genome_labels height_inches width_inches xmax ymax left_margin)) {
		die "Read_config: $_ not defined. Check or remake\n" if(!defined $$synima_data{$_});
	}
	die "Read_config: plot_synteny should equal c (contig) or g (gene)\n" if(($$synima_data{'plot_synteny'} ne 'c') && ($$synima_data{'plot_synteny'} ne 'g'));
	return 1;
}

sub make_array_from_config_line {
	my ($config_file, $key) = @_;
	my @array = split /,/, $$config_file{$key};
	$array[0] =~ s/^$key=//;
	return \@array;
}

sub make_hash_from_config_line {
	#my ($config_file, $key) = @_;
	my $config_line = $_[0];
	my %info;
	#my @array = split /=/, $$config_file{$key};
	my @array = split /=/, $config_line;
	die "make_hash_from_config_line called with wrong value: $config_line\n" if((!defined $array[1]) || (!defined $array[2]));
	$info{$array[1]} = $array[2];
	return \%info;
}

sub check_input_file_make_config_entry {
	my ($file, $type) = @_;
	my $cwd = getcwd;
	my $entry = '';
	
	if($file ne '') {
		if (! -e $file) { warn "Cannot open (skipping) file $file : $!\n"; } 
		else {
			if(-e "$cwd/$file") { $entry = "$type=$cwd/$file\n"; }
			else { $entry = "$type=$file\n"; }
		}
	}
	return $entry;
}

sub gather_genome_lengths_for_config {
	my ($input, $genomes_suffix, $verbose) = @_;

	# Save genome names
	my $genome_names_hash = daglines::save_genome_names_hash_from_aligncoords_spans($input);

	# Save the genome lengths from FASTA files
	my $genome_names;
	my @genome_names_array;
	my $genome_lengths;
	my $genome_length_max = 0;
	my $genome_contig_lengths;
	my %genome_contig_order;
	foreach my $genome(keys %{$genome_names_hash}) {
		my $fasta_file = "$genome/$genome.$genomes_suffix";
		die "Can't find the genome sequence ($fasta_file) to get full chromosome lengths: $!\n" if (! -e $fasta_file);

		# Save genome names array
		$genome_names .= "$genome,";
		push @genome_names_array, $genome;

		# Save genome -> lengths
		my $length = fastafile::fasta_to_total_seq_length($fasta_file);
		$genome_lengths .= "genome_lengths_$genome=$length\n";

		# Max length
		if($length > $genome_length_max) { $genome_length_max = $length; }

		# Save genome -> id -> lengths
		my $fasta_id_to_seq = fastafile::fasta_id_to_seq_length_hash($fasta_file);
		foreach my $id(sort keys %{$fasta_id_to_seq}) {
			my $length = $$fasta_id_to_seq{$id};
			$genome_contig_lengths .= "genome_contig_lengths_$genome=$id=$length\n";
		}

		# Save fasta contig order
		my ($order) = fastafile::fasta_id_to_order_array($fasta_file);
		$genome_contig_order{$genome} = $order;
	}
	$genome_names =~ s/,$//;
	$genome_lengths =~ s/\n$//;

	# Calculate number of genomes and half way point through genomes
	my $num_genomes = scalar(keys(%{$genome_names_hash}));
	if($verbose eq 'y') { warn "Found syntenic information for $num_genomes genomes\n"; }
	my $half_way_point_through_genomes = (($num_genomes / 2) -1);

	# Find contig order using 1st genome as guide
	my $genome_contig_order_string;
	for(my $i=0; $i<($num_genomes -1); $i++) {
		my $name1 = $genome_names_array[$i];
		my $name2 = $genome_names_array[($i + 1)];
		if($i eq 0) {
			$genome_contig_order_string = "contig_order_$name1=";
			foreach(@{$genome_contig_order{$name1}}) { $genome_contig_order_string .= "$_,"; }
			$genome_contig_order_string =~ s/,$//;
			$genome_contig_order_string .= "\n";
		}
		if($verbose eq 'y') { warn "find order of $name2 from $name1\n"; }
		my ($order2, $order2_array) = daglines::find_contig_order_from_aligncoords_spans_and_fastas($genome_contig_order{$name1}, $genome_contig_order{$name2}, $name1, $name2, $input, $verbose);
		$genome_contig_order{$name2} = $order2_array;
		$genome_contig_order_string .= "$order2\n";
	}
	$genome_contig_order_string =~ s/\n$//;

	return ($genome_names, $genome_lengths, $genome_contig_lengths, $num_genomes, $half_way_point_through_genomes, $genome_length_max, $genome_contig_order_string);
}

sub print_chromosome_name {
	my ($chr_name, $chr_length, $chr_name_pos_x, $chr_name_pos_y, $genome, $config_data, $ofh) = @_;
	return if($$config_data{'contig_min_length_for_label'} > $chr_length);

	# Replace terms
	#warn "Printing name for $chromosome_name\n";
	my $contig_name_replace_array = &make_array_from_config_line($config_data, 'contig_name_replace');
	my $chr_name_new = $chr_name;
	foreach(@{$contig_name_replace_array}) { 
		my @parts = split /:/, $_;
		$chr_name_new =~ s/^$parts[0]/$parts[1]/;
	}

	# Reverse complimented?
	my $chr_rev = daglines::genome_to_chromosome_to_length_hash_subset($config_data);
	if(defined $$chr_rev{$genome}{$chr_name}) { $chr_name_new .= "-"; }
	else { $chr_name_new .= "+"; }

	# Print and return
	print $ofh "text($chr_name_pos_x, $chr_name_pos_y, \"$chr_name_new\", cex = 1.0, srt=90)\n";
	return;
}

sub plot_genes {
	my ($up, $down, $cumulative_length, $transparency, $rgb, $lwd, $lty, $points_or_polygon, $genes, $ofh) = @_;
	GENES: foreach my $gene_ids(keys %{$genes}) {
		GENES: foreach my $gene_start(keys %{$$genes{$gene_ids}}) {
			my $gene_stop = $$genes{$gene_ids}{$gene_start};
			my $gene_start_cumulative = ($gene_start + $cumulative_length);
			my $gene_stop_cumulative = ($gene_stop + $cumulative_length);
			if($points_or_polygon eq 'polygon') { 
				print $ofh "x <- c($gene_start_cumulative, $gene_start_cumulative, $gene_stop_cumulative, $gene_stop_cumulative, $gene_start_cumulative)\n";
				print $ofh "y <- c($down, $up, $up, $down, $down)\n";
				print $ofh "polygon(x, y, lty = $lty, lwd = $lwd, col=adjustcolor(\"$rgb\", alpha.f=1), border=\"NA\")\n"; 
			}
			elsif($points_or_polygon eq 'points') { print $ofh "points($gene_start_cumulative, y = $up, pch=19, col=adjustcolor(\"$rgb\", alpha.f=1), cex=0.8)\n"; } 
			else { die "$points_or_polygon neither points or polygon\n"; }
		}
	}
	return 1;
}

sub print_chromosome_synteny {
	my ($ss1, $ss2, $y1, $rgb, $ofh) = @_;
	my $y2 = ($y1+1);
	print $ofh "x <- c($$ss1[0], $$ss2[0], $$ss2[1], $$ss1[1], $$ss1[0])\n";
	print $ofh "y <- c($y1, $y2, $y2, $y1, $y1)\n";
	print $ofh "polygon(x, y, lty = 2, lwd = 2, col=adjustcolor(\"$rgb\", alpha.f=0.5), border=NA)\n";
	return 1;
}

sub print_gene_synteny {
	my ($ss1, $ss2, $i, $colour, $ofh) = @_;
	my $syn_up = ($i+1);
	my $transparency = 0.5;
	my ($lwd, $lty) = (1, 1);
	print $ofh "x <- c($$ss1[0], $$ss2[0], $$ss2[1], $$ss1[1], $$ss1[0])\n";
	print $ofh "y <- c($i, $syn_up, $syn_up, $i, $i)\n";
	print $ofh "polygon(x, y, lty = $lty, lwd = $lwd, col=adjustcolor(\"$colour\", alpha.f=$transparency), border=NA)\n";
	return 1;
}

sub find_y_pos_for_genes {
	my ($gene_up, $gene_down, $half_way_point, $up1, $down1, $up2, $down2) = @_;
	if($gene_up <= $half_way_point) {
		$gene_up += $up1;
		$gene_down = $gene_up-$down1;
	} else {
		$gene_up -= $up2;
		$gene_down = $gene_up-$down2;
	}
	return ($gene_up, $gene_down);
}

sub determine_gene_synteny_colour_and_positions {
	my ($synima_data, $genome_synteny, $gene_start_stop1, $gene_start_stop2, $cumulative_length, $chr2, $chr_order2, $genome_contig_lengths2, $name1, $name2, $chromosomes) = @_;

	# Cumulative positions
	my ($gene_ss1, $gene_ss2) = &find_cumulative_positions($gene_start_stop1, $gene_start_stop2, $cumulative_length, $chr2, $chr_order2, $genome_contig_lengths2);

	# Colour
	my $colour = $$synima_data{'aligncoords_spans_colour1'};
	for(my $i=2; $i<=3; $i++) {
		my $filename = "aligncoords_spans_location$i";
		my $colour_id = "aligncoords_spans_colour$i";
		my $presence = daglines::check_hash_keys_from_aligncoords($$genome_synteny{$filename}, 4, $name1, $name2, $chromosomes, $chr2); 
		if($presence eq 1) {
			foreach my $start_stop1(keys %{$$genome_synteny{$filename}{$name1}{$name2}{$chromosomes}{$chr2}}) {
				my $start_stop2 = $$genome_synteny{$filename}{$name1}{$name2}{$chromosomes}{$chr2}{$start_stop1};
				my ($ss1, $ss2) = &find_cumulative_positions($start_stop1, $start_stop2, $cumulative_length, $chr2, $chr_order2, $genome_contig_lengths2);
				#warn "CHECK IF MY CURRENT GENE $$gene_ss1[0] - $$gene_ss1[1] and $$gene_ss2[0] - $$gene_ss2[1] is inside $$ss1[0] - $$ss1[1] and $$ss2[0] - $$ss2[1]\n";
				if((($$gene_ss1[0] >= $$ss1[0]) && ($$gene_ss1[0] <= $$ss1[1])) && (($$gene_ss2[0] >= $$ss2[0]) && ($$gene_ss2[0] <= $$ss2[1]))) { $colour = $$synima_data{$colour_id}; }
			}
		}
	}
	return($colour, $gene_ss1, $gene_ss2);
}

sub find_cumulative_positions {
	my($sns1, $sns2, $sns1_cumulative, $sns2_chr, $sns2_chr_order, $genome2_lengths) = @_;
	my @ss1 = split /-/, $sns1;
	my @ss2 = split /-/, $sns2;
	
	# cumulative position
	foreach(@ss1) { $_ += $sns1_cumulative; }
	
	# Find the cumulative position of my sytenic genome above
	my $cumulative_length2 = 0;
	POSITIONGENOMETWO: foreach my $chs(@{$sns2_chr_order}) {
		last POSITIONGENOMETWO if($chs eq $sns2_chr);
		die "Problem encountered in find_cumulative_positions for $chs. Re-run\n" if(!defined $$genome2_lengths{$chs});
		$cumulative_length2 += $$genome2_lengths{$chs};
	}
	foreach(@ss2) { $_ += $cumulative_length2; }
	#warn "cumulatives = @ss1 - @ss2\n";
	return(\@ss1, \@ss2);
}

sub process_cmd {
	my ($cmd) = @_;
	warn "CMD: $cmd\n";
	my $ret = system($cmd);
	die "Error, cmd $cmd died with return $ret\n" if($ret);
	return 1;
}

1;
