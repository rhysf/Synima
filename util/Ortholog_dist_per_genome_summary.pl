#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules";
use read_Tab;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <cluster_dist_per_genome.txt> > output_num_genes.tab\n
Optional: -c\tclusters_and_unique.wAnnots
          -l\tSeperate file with info to join []
	  -d\tIf opt_l, then from which column do i look for id? [0]
          -p\tPrinting options (n=none, c=look up names in opt_c and print those lines, l=get gene names from opt_c, then look up in opt_l and print those lines) [n]
	  -z\tIf opt_p = c, then restrict to 1 description - E.g. missing_in_Fo5176_454 []\n
Notes: Output files go to -c (wAnnots) directory\n";
our($opt_r, $opt_c, $opt_l, $opt_d, $opt_p, $opt_z);
getopt('rcldpz');
die $usage unless ($opt_r);
foreach($opt_r) { die "file $_ is not readable: $!\n" if (! -e $_); }
if(!defined $opt_p) { $opt_p = 'n'; }
if(!defined $opt_d) { $opt_d = 0; }
if($opt_p eq 'l') { 
	die "file -l needs to be specified for -p l\n" if(!defined $opt_l);
	die "file -l is not readable: $!\n" if(! -e $opt_l);
}

# Save info from wAnnots file and/or additional file (column d -> line)
my ($wAnnots, $file_l);
if(($opt_p eq 'c') || ($opt_p eq 'l')) { $wAnnots = &save_wAnnots($opt_c); }
if($opt_p eq 'l') { $file_l = &save_additional_file($opt_l, $opt_d); }

# How many isolates in cluster dist
my ($number_of_isolates, $isolate_names) = &save_isolate_names_and_count_from_cluster_dist($opt_r);

my %num_genes;
my %num_clusters;
my %clusters;
my %missing_counts;

# Go through ortholog file
open my $fh, '<', $opt_r or die "Cannot open $opt_r : $!\n";
GENEORTHOLOGS: while(my $line=<$fh>) {
	chomp $line;
	next GENEORTHOLOGS if($line =~ m/^\#/);
	my @bits = split /\t/, $line;

	# Grab only the gene counts
	my @genes_in_cluster = @bits[ 2 .. $#bits ];
	my $counts_of_genes_in_cluster = join ("\t", @genes_in_cluster);
	my $sum_of_genes_in_cluster = 0;
	for (@genes_in_cluster) { $sum_of_genes_in_cluster += $_; }
	my $number_of_isolates_minus_1 = ($number_of_isolates -1);
	my $desc;

	### ORDER IS IMPORTANT
	# missing in something
	if($counts_of_genes_in_cluster =~ m/^0\t|\t0\t|\t0$/)  {
		
		#warn "$counts_of_genes_in_cluster missing in something\n";
		$desc = 'missing_in_';
		my $tmp_missing_count = 0;
		my $tmp_isolate_name;
		for(my $i=2; $i<scalar(@bits); $i++) {
			if($bits[$i] eq 0) { 
				$desc .= $$isolate_names{$i} . "_and_"; 
				$missing_counts{$$isolate_names{$i}}++;
				$tmp_missing_count++;
			}
		}
		$desc =~ s/\_and\_$//;

		# new 2023... for large numbers of comparisons, this description gets unweildy.
		# it would be good to see if instead of 'missing in', i could include some 'unique_in'
		# especially if it's lost in everything except 1 genome
		if($tmp_missing_count eq $number_of_isolates_minus_1) {
			#warn "unique in 1: $line\n";
			for(my $i=2; $i<scalar(@bits); $i++) {
				if($bits[$i] eq 1) { 
					$desc = 'unique_in_' . $$isolate_names{$i};
				}
			}
		}
	}
	# Not missing in anything, paralog in something
	elsif($counts_of_genes_in_cluster =~ m/[23456789]|10|11/) {
		#warn "$counts_of_genes_in_cluster paralog in something\n";
		$desc = 'paralog_in_';
		# single copy in anything?
		if($counts_of_genes_in_cluster !~ m/^1\t|\t1\t|\t1$/) { $desc .= 'all'; } 
		else {
			for(my $i=2; $i<scalar(@bits); $i++) {
				if($bits[$i] > 1) { $desc .= $$isolate_names{$i} . "_and_"; }
			}
			$desc =~ s/\_and\_$//;
		}
	}
	# 1:1 in all
	elsif($counts_of_genes_in_cluster =~ m/^(1\t){$number_of_isolates_minus_1}1$/) {
		#warn "$counts_of_genes_in_cluster is 1:1\n";
		$desc = 'orthologs';
	}
	else {
		warn "Have not accounted for this line: $line\n";
	}
	
	# Summary
	$num_genes{$desc}+=$sum_of_genes_in_cluster; 
	$num_clusters{$desc}++;
	$clusters{$desc}{$bits[0]} = 1;
}

# Additional summary of counts
warn "Summary 1:\nIsolate\tClusters_missing\n";
foreach my $genome(sort keys %missing_counts) {
	warn "$genome\t$missing_counts{$genome}\n";
}

## Summary of counts
print "Type\tClusters\tGenes\n";
foreach my $type(sort keys %num_genes) {
	my $tally_clust = $num_clusters{$type};
	my $tally_genes = $num_genes{$type};
	print "$type\t$tally_clust\t$tally_genes\n";
}

# l=get gene names from opt_c, then look up in opt_l and print those lines
if(($opt_p eq 'c') || ($opt_p eq 'l')) {
	warn "finding in wAnnots file...\n";
	DESC: foreach my $desc(keys %clusters) {
		## Specify output
		my $out;
		if($opt_p eq 'c') { $out = "$opt_c-" . $desc; }
		else { $out = "$opt_l-" . $desc; }

		if($opt_z) {
			next DESC if($opt_z ne $desc);
		}

		my $count = 0;
		open my $ofh, '>', $out or die "Cannot open $out: $!\n";
		foreach my $cluster_number(keys %{$clusters{$desc}}) {
			die "haven't got anything saved from wAnnots file for $cluster_number\n" if(!defined $$wAnnots{$cluster_number});
			foreach my $long_gene_id(keys %{$$wAnnots{$cluster_number}}) {

				# Check i've got the right IDs
				if($count eq 0) {
					$count = 1;
					warn "The kind of IDs i'm saving = $long_gene_id . If this doesn't look right... edit script or parameters to specify correct gene id\n";
				}

				# either print or save long id for whatever is in opt_l
				if($opt_p eq 'l') { 
					if(defined $$file_l{$long_gene_id}) { 
						print $ofh $$file_l{$long_gene_id} . "\n"; 
					}
				} else {
					my $wAnnots_line = $$wAnnots{$cluster_number}{$long_gene_id};
					#die "check $cluster_number -> $long_gene_id -> $wAnnots_line\n";
					print $ofh "$wAnnots_line\n";
				}
			}
		}
		close $ofh;
	}
}

sub save_additional_file {
	my ($file, $column) = @_;
	my %info;
	open my $fh, '<', $file or die "Can't open $file : $!\n";
	FILE: while(my $line = <$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		$info{$bits[$column]} = $line;
	}
	close $fh;
	return \%info;
}

sub save_wAnnots {
	my ($file) = $_[0];
	my %cluster_to_long_id;
	open my $fh, '<', $file or die "Can't open $file : $!\n";
	GENEINFO2: while(my $line = <$fh>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($cluster,$genome,$long_id,$parent_id,$short_id,$annotation,$contig,$start,$stop,$strand,$fams) = @bits;
		next GENEINFO2 if($line =~ m/^\n/);
		next GENEINFO2 if(!defined $bits[0]);
		$cluster_to_long_id{$cluster}{$parent_id} = $line;
	}
	close $fh;
	return \%cluster_to_long_id;
}

sub save_isolate_names_and_count_from_cluster_dist {
	my $file = $_[0];
	my $number_of_isolates = 0;
	my %isolate_names;
	open my $fh, '<', $file or die "Cannot open $file : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		if($line =~ m/^#cluster_id/) {
			my @bits = split /\t/, $line;
			$number_of_isolates = (scalar(@bits) - 2);
			warn "$number_of_isolates number of isolates found\n";
			for(my $i=2; $i<scalar(@bits); $i++) { $isolate_names{$i} = $bits[$i]; warn "$bits[$i]\n"; }
		}
		elsif($line =~ m/^#/) { next; }
		else { die "$file does not have the expected cluster.dist header : $line\n"; }
		last;
	}
	close $fh;
	return ($number_of_isolates, \%isolate_names);
}
