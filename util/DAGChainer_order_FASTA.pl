#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/perl_modules";
use MutationTools::read_FASTA;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: $0 -a <current_order_isolate1> -b <name_isolate2> -s <aligncoords.spans> > new_order_isolate2\n
Optional: -f fasta for isolate2 (will include everything that has no synteny at the end of the new_order_isolate2\n";
our($opt_a, $opt_b, $opt_s, $opt_f);
getopt('absf');
die $usage unless ($opt_a && $opt_b && $opt_s);

# save current order
my ($order1, $isolate1) = &save_order($opt_a);
my $fasta_order2;
if(defined $opt_f) { 
	my ($sequences, $descriptions, $order) = fastafile::fasta_id_to_seq_hash($opt_f);
	$fasta_order2 = $order;
}

# find matches
my @new_order2;
my %contigs_found;
ORDER: foreach my $chr_order(@{$order1}) {
	my $search_for = "grep \$'$isolate1;$chr_order\\t' $opt_s | grep $opt_b";
	warn "Looking for syntenic regions between $isolate1 ($chr_order) and $opt_b... using $search_for\n";
	my $lines = `$search_for`;
	die "found nothing for: grep \$'$isolate1;$chr_order\t' $opt_s | grep $opt_b.... investigate.\n" if(!defined $lines);
	#next ORDER if(!defined $lines);
	my @line = split /\n/, $lines;
	my (%new_order_info, %align_lengths);
	foreach my $chr_line(@line) {
		my @bits = split /\t/, $chr_line;
		my ($is1, $ss1, $length1, $is2, $ss2, $length2) = @bits;
		my @isolate_and_chr1 = split /;/, $is1;
		my @isolate_and_chr2 = split /;/, $is2;
		my @start_stop1 = split /-/, $ss1;
		my @start_stop2 = split /-/, $ss2;

		# find order
		# I want order from isolate1, and then to save isolate2 chromosome names. 
		if($isolate_and_chr1[0] eq $isolate1) {	$new_order_info{$start_stop1[0]}{$start_stop1[1]} = $isolate_and_chr2[1]; }
		if($isolate_and_chr2[0] eq $isolate1) {	$new_order_info{$start_stop2[0]}{$start_stop2[1]} = $isolate_and_chr1[1]; }
		
		# Length (can be accumulative here) should be used to decide if multiple hits are found
		$align_lengths{$isolate_and_chr1[1]} += $length1;
		$align_lengths{$isolate_and_chr2[1]} += $length1;
	}

	# sort them and save order. Still has redundancies
	foreach my $start(sort { $a <=> $b } keys %new_order_info) {
		foreach my $stop(keys %{$new_order_info{$start}}) {
			my $contig = $new_order_info{$start}{$stop};

			# check if it is already there?
			if(defined $contigs_found{$contig}) {
			
				# if the last time it was found, it was smaller - then replace with this one
				if($contigs_found{$contig} < $align_lengths{$contig}) {
					# delete it from array
					for(my $i=0; $i<scalar(@new_order2); $i++) {
						if($new_order2[$i] eq $contig) { splice @new_order2, $i, 1; }
					}
					#warn "$contig moved to front due to $contigs_found{$contig} < $align_lengths{$contig}\n";

					# replace with new value
					push @new_order2, $contig;
					$contigs_found{$contig} = $align_lengths{$contig};
				}
			} else {
				my $length = $align_lengths{$contig};
				$contigs_found{$contig} = $align_lengths{$contig};
				push @new_order2, $contig;
				#warn "putative order (based on $chr_order): $contig\n";
			}
		}
	}
}

# print new file
my %found;
print "#$opt_b\n";
foreach(@new_order2) { print "$_\n"; $found{$_}=1; }

# Add missing chromosomes/contigs etc.
ORDER: foreach my $chr_order(@{$fasta_order2}) {
	if(!defined $found{$chr_order}) { print "$chr_order\n"; }
}

sub save_order {
	my $file = $_[0];
	open IN, "<$file" or die "Cannot read $file: $!\n";
	my $genome_found;
	my @order;
	LINES: while(my $line=<IN>) {
		chomp $line;
		if($line =~ m/#/) {
			$line =~ s/^#//;
			$genome_found = $line;
			next LINES;
		} 
		if(!defined $genome_found) { die "$file does not start with #isolate_name\n"; }
		push @order, $line;
	}
	return (\@order, $genome_found);
}
