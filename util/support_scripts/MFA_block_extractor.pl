#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../../modules/";
#use Fasta_reader;
use read_FASTA;
use Data::Dumper;
use List::Util qw (min max);

### r.farrer@exeter.ac.uk (edited from Chris D. 8/3/12)

# Opening statement
my $usage = "Usage: perl $0 <file.msa> [MIN_BLOCK_LENGTH=20] [BLOCK_EDGE_EXCLUDE=0]\n";
die $usage unless(@ARGV >= 1);
my $malign_file = $ARGV[0];
my $MIN_BLOCK_LENGTH = $ARGV[1] || 20;
my $BLOCK_EDGE_EXCLUDE = $ARGV[2] || 0; ## changed default to 0


my @entries;
#my $fasta_reader = new Fasta_reader($malign_file);
my $fasta = fastafile::fasta_to_struct($malign_file);

foreach my $id(sort keys %{$$fasta{'seq'}}) {
	my $seq = $$fasta{'seq'}{$id};
	my $align_struct = &get_alignment_struct($id, $seq);
	push (@entries, $align_struct);
}
	
# get bounds for all extraction
my @lends;
my @rends;
foreach my $entry (@entries) {
	my $lend = $entry->{lend};
	my $rend = $entry->{rend};
		
	push (@lends, $lend);
	push (@rends, $rend);
}

my $lend_bound = max(@lends);
my $rend_bound = min(@rends);

die "Error, no contiguous alignment block is possible.  Bound rend: $rend_bound < bound lend: $lend_bound" if ($rend_bound < $lend_bound);


my $block_start = undef;
for (my $i = $lend_bound; $i <= $rend_bound; $i++) {
		
	my $block_term_flag = 0;
	if (&all_columns_populated(\@entries, $i)) {
		# start or extend block
		if (! defined $block_start) {
			$block_start = $i;
		}
	}
	else {
		$block_term_flag = 1;
	}

	if ($block_term_flag || $i == $rend_bound) {

		if (defined $block_start &&  $i - $block_start >= $MIN_BLOCK_LENGTH) {
			my $block_end = ($i == $rend_bound && ! $block_term_flag) ? $i : $i - 1; # include last column if fully populated
			&append_block(\@entries, $block_start, $block_end);
		}
		$block_start = undef; # reset for next block
	}
}

# print blocks
foreach my $entry (@entries) {
	my $acc = $entry->{acc};
	my $blocks = $entry->{blocks};

	$blocks =~ s/(\S{60})/$1\n/g;
	chomp $blocks;
		
	print ">$acc\n$blocks\n";
}

sub all_columns_populated {
	my ($entries_aref, $column) = @_;

	foreach my $entry (@$entries_aref) {
		
		my $char = $entry->{chars}->[$column];
		if ($char !~ /\w/) {
			return(0); # not populated (gap instead)
		}
	}

	## all entries populated
	return(1);
}

sub append_block {
	my ($entries_aref, $block_start, $block_end) = @_;
	
	$block_start += $BLOCK_EDGE_EXCLUDE;
	$block_end -= $BLOCK_EDGE_EXCLUDE;
	
	foreach my $entry (@$entries_aref) {
		for (my $i = $block_start; $i <= $block_end; $i++) {
			my $char = $entry->{chars}->[$i];
			$entry->{blocks} .= $char;
			
			## add spacer if last position of block. Changed "." to "-" to prevent RAxML char errors, Paul G. 8/3/2012
			#if ($i == $block_end) { ## Removed block separator, Chris D. 8/3/2012
			#	$entry->{blocks} .= "-";
			#}
			

		}
	}
	return;
}

sub get_alignment_struct {
	my ($acc, $alignment) = @_;

	my @chars = split(//, $alignment);
	my $lend = &find_left_alignchar(@chars);
	my $rend = &find_right_alignchar(@chars);
	my $num_residues = scalar (grep { /[A-Za-z]/ } @chars);
	
	my $struct = { 
		chars => \@chars,
		acc => $acc,
		lend => $lend,
		rend => $rend,
		num_residues => $num_residues,
		blocks => "",
	};
	return($struct);
		
}

sub find_left_alignchar {
	my @chars = @_;
	for (my $i = 0; $i <= $#chars; $i++) {
		if ($chars[$i] =~ /[A-Za-z]/) {
			return($i);
		}
	}
	die "Error, no aligned char found from lend: @chars";
}

sub find_right_alignchar {
	my @chars = @_;
	for (my $i = $#chars; $i >= 0; $i--) {
		if ($chars[$i] =~ /[A-Za-z]/) {
			return($i);
		}
	}
	die "Error, no aligned char found from rend: @chars";
}
