#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use Cwd;
use File::Basename;
use Data::Dumper;
use FindBin qw($Bin);
use lib "$Bin/modules";
use Synima;
#use read_FASTA;
use read_DAG;
use read_Tab;
use write_Rscript;
my $cwd = getcwd;

### rfarrer@broadinstitute.org

# Opening commands 
my $usage = "Usage:     perl $0 -c <config.txt> or -a <aligncoords> -b <aligncoords.spans>\n
Commands:  -c\tConfig.txt [$cwd/SynIma-output/config.txt]
           -a\tAligncoords
           -b\tAligncoords.spans\n
Optional:  -e\tGenome FASTA filename extension (e.g. $cwd/genome1/genome1.genome.fa etc.) [genome.fa]
           -t\tAligncoords.spans 2
	   -u\tAligncoords.spans 3
	   -k\tGene IDs 1 (1 per line)
	   -l\tGene IDs 2 (1 per line)
	   -o\tGene IDs 3 (1 per line)
	   -r\tRun full program (y) or just create config (n) [y]
	   -v\tVerbose output (y/n) [n]\n
Plot Opts: -i\tWidth of figure in pixels [1100]
	   -j\tHeight of figure in pixels (num of genomes * 100)
           -g\tFill in chromosome/contig synteny (c) or gene synteny (g) [c]
	   -z\tPlot individual genes (y/n) [n]
           -x\tOrder of genomes from bottom to top seperated by comma
	   -n\tGenome labels from bottom to top seperated by comma
	   -w\tnumber of lines for left hand margin [12]\n
Notes:     Config.txt will be made automatically if not present, and read automatically if it is.
           Config.txt specifies order of genomes, chromosomes, colours, and other plot options. 
	   Config.txt can be manually edited after creation.
	   Default genome labels will be as they appear in aligncoords
	   Order of genomes must have names as they appear in aligncoords
	   Aligncoords.spans and Gene ID files will be highlighted according to the config\n";
our($opt_a, $opt_b, $opt_c, $opt_d, $opt_e, $opt_f, $opt_g, $opt_h, $opt_i, $opt_j, $opt_k, $opt_l, $opt_m, $opt_n, $opt_o, $opt_p, $opt_q, $opt_r, $opt_s, $opt_t, $opt_u, $opt_v, $opt_w, $opt_x, $opt_y, $opt_z);
getopt('abcdefghijklmnopqrstuvwxyz');
die $usage unless (($opt_c) || ($opt_a && $opt_b));
foreach($opt_a, $opt_b, $opt_k, $opt_l, $opt_n, $opt_o, $opt_t, $opt_u) { if(!defined $_) { $_ = ''; } }
if(!defined $opt_c) { $opt_c = "$cwd/SynIma-output/config.txt"; }
if(!defined $opt_e) { $opt_e = 'genome.fa'; }
if(!defined $opt_g) { $opt_g = 'c'; }
if(!defined $opt_i) { $opt_i = 1100; }
if(!defined $opt_j) { $opt_j = 'n'; }
if(!defined $opt_r) { $opt_r = 'y'; }
if(!defined $opt_v) { $opt_v = 'n'; }
if(!defined $opt_w) { $opt_w = 12; }
if(!defined $opt_x) { $opt_x = 'n'; }
if(!defined $opt_z) { $opt_z = 'n'; }
die "Option g must equal c or g : $opt_g\n" if(($opt_g ne 'c') && ($opt_g ne 'g'));
warn "Running Synima...\n";

# Read Config.txt or make one
my $synima_data;
if(-e $opt_c) { $synima_data = synima::read_config($opt_c); } 
else { $synima_data = synima::write_config($opt_a, $opt_b, $opt_c, $opt_e, $opt_g, $opt_i, $opt_j, $opt_k, $opt_l, $opt_n, $opt_o, $opt_r, $opt_t, $opt_u, $opt_v, $opt_w, $opt_x, $opt_z); }
die "Finished reading/writing config file. Run specified to no\n" if($opt_r eq 'n');

# Save gene synteny
my ($gene_synteny, $genes) = daglines::save_aligncoords_and_reverse_wrapper_from_config($synima_data);

# Save genome synteny and gene subsets (specified as lists)
my $genome_synteny;
my $genes_subsets;
for(my $i=1; $i<=5; $i++) {
	my $file_name1 = "aligncoords_spans_location$i";
	my $file_name2 = "gene_list$i";
	if($$synima_data{$file_name1}) { $$genome_synteny{$file_name1} = daglines::save_aligncoords_spans_and_reverse_wrapper_from_config($synima_data, $file_name1); }
	if($$synima_data{$file_name2}) { $$genes_subsets{$file_name2} = daglines::genes_subset_from_gene_list_from_config($synima_data, $file_name2, $genes); }
}

# Prepare variables for plot
my $order = synima::make_array_from_config_line($synima_data, 'genome_order');
my $ylab = synima::make_array_from_config_line($synima_data, 'genome_labels');
die "Different number of genomes in genome_order and genome_labels. Remake or edit $opt_c\n" if(scalar(@{$order}) ne scalar(@{$ylab}));
#die "Different number of genomes in aligncoords_spans_location1 ($num_genomes_saved_in_aligncoords1) than specified in order or ylab\n" if($num_genomes_saved_in_aligncoords1 ne (scalar(@{$order})));
my ($running_genome_length);
foreach my $genome_id(@{$order}) { $$running_genome_length{$genome_id} = 0; }

# Begin R file
if($opt_v eq 'y') { warn "Creating Rscript $$synima_data{'output_rscript'}...\n"; }
open my $ofh, '>', $$synima_data{'output_rscript'} or die "Cannot open" . $$synima_data{'output_rscript'} . ": $!\n";
print $ofh "#! /usr/bin/Rscript --vanilla\n";
print $ofh "suppressMessages(expr = \"SynIma: Generating pdf. May take some time...\")\n";
print $ofh "require(graphics)\n";
print $ofh "pdf(\"$$synima_data{'output_pdf'}\", height=$$synima_data{'height_inches'}, width=$$synima_data{'width_inches'})\n";
print $ofh "par(mfrow=c(1,1), mar=c(0, $$synima_data{'left_margin'}, 0, 0), oma=c(4, 0, 4, 0), las=1, cex.axis=1.2)\n";
print $ofh "plot(0,0, ylab=\"\", yaxt=\"n\", xlim=c(0, $$synima_data{'xmax'}), xaxt = \"n\", ylim=c(-0.1,$$synima_data{'ymax'}), cex.lab=2, xaxs = \"r\")\n";
print $ofh "box(which = \"plot\", lty = \"solid\", lwd=3, col=\"white\")\n";

# Begin R file plot
for(my $i=0; $i<$$synima_data{'genome_number_of'}; $i++) {

	# xlab and xaxis on the bottom
	if($i eq ($$synima_data{'genome_number_of'} - 1)) { rlines::print_xlab_and_xaxis_for_windows($$synima_data{'xmax'}, $ofh); }

	# Names and labels for y-axis
	my $name1 = $$order[$i];
	my $name2 = $$order[($i+1)];
	if(!defined $name2) { $name2 = $$order[0]; }
	if($opt_v eq 'y') { warn "Plotting synteny for $name1 - $name2...\n"; }
	print $ofh "axis(2, at=$i, lwd.ticks=0, lab=c(\"$$ylab[$i]\"))\n";

	# Chromosome/contig lengths and order
	foreach my $name($name1, $name2) { 
		die "Bad genome_contig_lengths_$name. Remake or edit $opt_c\n" if(! $$synima_data{"genome_contig_lengths_$name"}); 
		die "Bad contig_order_$name. Remake or edit $opt_c\n" if(! $$synima_data{"contig_order_$name"}); 
	}
	my $genome_contig_lengths1 = $$synima_data{"genome_contig_lengths_$name1"};
	my $genome_contig_lengths2 = $$synima_data{"genome_contig_lengths_$name2"};
	my $chr_order1 = synima::make_array_from_config_line($synima_data, "contig_order_$name1");
	my $chr_order2 = synima::make_array_from_config_line($synima_data, "contig_order_$name2");

	# Positions on Chromosomes/contigs
	CHROMOSOMES: foreach my $chromosomes(@{$chr_order1}) {
		if($opt_v eq 'y') { warn "Defining positions on contig $chromosomes\n";	}

		# Contig positions
		die "genome_contig_lengths_$name1=$chromosomes not found in config. Remake or edit $opt_c\n" if(!defined $$genome_contig_lengths1{$chromosomes});
		die "genome_contig_lengths_$name1=$chromosomes corrupted. Remake or edit $opt_c\n" if($$genome_contig_lengths1{$chromosomes} !~ m/\d+/);
		my $length = $$genome_contig_lengths1{$chromosomes};
		my $half_contig_length = ($length / 2);
		my $cumulative_length = $$running_genome_length{$name1};
		my $length_plus_cumulative = ($length + $cumulative_length);
		my $contig_half_way = ($length_plus_cumulative - $half_contig_length);

		# Contig lines and name (text above)
		my $text_y = ($i + 0.3);
		my $up = ($i+0.03);
		my $down = ($i-0.03);
		print $ofh "segments(x0 = $cumulative_length, y0 = $i, x1 = $length_plus_cumulative, y1 = $i, col = \"black\", lty=1)\n";
		print $ofh "segments(x0 = $cumulative_length, y0 = $down, x1 = $cumulative_length, y1 = $up, col = \"black\", lty=1, lwd=2)\n";
		print $ofh "segments(x0 = $length_plus_cumulative, y0 = $down, x1 = $length_plus_cumulative, y1 = $up, col = \"black\", lty=2, lwd=3)\n";
		synima::print_chromosome_name($chromosomes, $length, $contig_half_way, $text_y, $name1, $synima_data, $ofh);

		# Contig synteny (only for genomes before the final as no more synteny afterwards)
		if(($$synima_data{'plot_synteny'} eq 'c') && ($i < ($$synima_data{'genome_number_of'} - 1))) {
			foreach my $chr2(keys %{$$genome_synteny{'aligncoords_spans_location1'}{$name1}{$name2}{$chromosomes}}) {
				foreach my $start_stop1(keys %{$$genome_synteny{'aligncoords_spans_location1'}{$name1}{$name2}{$chromosomes}{$chr2}}) {
					my $start_stop2 = $$genome_synteny{'aligncoords_spans_location1'}{$name1}{$name2}{$chromosomes}{$chr2}{$start_stop1};	
					if($opt_v eq 'y') { warn "$name1 $chromosomes $start_stop1 - $name2 $chr2 $start_stop2\n"; }
					my ($ss1, $ss2) = synima::find_cumulative_positions($start_stop1, $start_stop2, $cumulative_length, $chr2, $chr_order2, $genome_contig_lengths2);

					# Determine colour and print
					my $colour = $$synima_data{'aligncoords_spans_colour1'};
					my $spans2_presence = daglines::check_hash_keys_from_aligncoords($$genome_synteny{'aligncoords_spans_location2'}, 5, $name1, $name2, $chromosomes, $chr2, $start_stop1); 
					my $spans3_presence = daglines::check_hash_keys_from_aligncoords($$genome_synteny{'aligncoords_spans_location3'}, 5, $name1, $name2, $chromosomes, $chr2, $start_stop1); 
					if($spans2_presence eq 1) { $colour = $$synima_data{'aligncoords_spans_colour2'}; }
					if($spans3_presence eq 1) { $colour = $$synima_data{'aligncoords_spans_colour3'}; }
					synima::print_chromosome_synteny($ss1, $ss2, $i, $colour, $ofh);
				}
			}
		}

		# Extend the chromosome running length
		$$running_genome_length{$name1} = $length_plus_cumulative;

		# Plot all genes?
		my ($gene_up, $gene_down) = synima::find_y_pos_for_genes($i, $i, $$synima_data{'genome_half_way_point'}, 0.03, 0.01, 0.02, 0.01);
		if($$synima_data{'plot_individual_genes'} eq 'y') {
			if($opt_v eq 'y') { warn "Plotting genes for $name1 $chromosomes...\n"; }
			synima::plot_genes($gene_up, $gene_down, $cumulative_length, 1, 'black', 2, 1, 'polygon', \%{$$genes{$name1}{$chromosomes}}, $ofh);
		}

		# Plot genes specified in lists
		GENELIST: for(my $i=1; $i<=5; $i++) {
			my $file_name = "gene_list$i";
			my $colour = "gene_list_colour$i";
			last GENELIST if(!defined $$synima_data{$colour});
			if($opt_v eq 'y') { warn "Plotting $file_name ($$synima_data{$colour})...\n"; }
			synima::plot_genes($gene_up, $gene_down, $cumulative_length, 1, $$synima_data{$colour}, 16,1, 'points', $$genes_subsets{$file_name}{$name1}{$chromosomes}, $ofh);
			($gene_up, $gene_down) = synima::find_y_pos_for_genes($gene_up, $gene_down, $$synima_data{'genome_half_way_point'}, 0.12, 0.03, 0.12, 0.03);
		}
		next CHROMOSOMES if($i >= ($$synima_data{'genome_number_of'} - 1));

		# Gene synteny wanted?
		next CHROMOSOMES if($$synima_data{'plot_synteny'} ne 'g');
		if($opt_v eq 'y') { warn "Plotting gene synteny...\n"; }
		foreach my $chr2(keys %{$$gene_synteny{$name1}{$name2}{$chromosomes}}) {
			if($opt_v eq 'y') { warn "Synteny between $name1 - $name2 on $chromosomes - $chr2\n"; }
			GSS1: foreach my $gene_start_stop1(keys %{$$gene_synteny{$name1}{$name2}{$chromosomes}{$chr2}}) {
				my $gene_start_stop2 = $$gene_synteny{$name1}{$name2}{$chromosomes}{$chr2}{$gene_start_stop1};
				if($opt_v eq 'y') { warn "Gene: $name1 $chromosomes $gene_start_stop1 - $name2 $chr2 $gene_start_stop2\n"; }

				# Determine colour and print
				my ($colour, $gene_ss1, $gene_ss2) = synima::determine_gene_synteny_colour_and_positions($synima_data, $genome_synteny, $gene_start_stop1, $gene_start_stop2, $cumulative_length, $chr2, $chr_order2, $genome_contig_lengths2, $name1, $name2, $chromosomes);
				synima::print_gene_synteny($gene_ss1, $gene_ss2, $i, $colour, $ofh);
			}
		}
	}
}

# Finish Rscript and run
rlines::close_device_print($ofh);
rlines::launch_R($$synima_data{'output_rscript'});
