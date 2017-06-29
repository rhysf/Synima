package roundnumber;
use strict;
use warnings;
use Math::Round;
use Data::Dumper;

### rfarrer@broadinstitue.org

sub round_number_to_nearest_number {
	my ($number, $round_to) = @_;

	# Nearest
	my $rounded_number = nearest($round_to, $number);
	return $rounded_number;
}

sub split_and_round_numbers_to_nearest_number {
	my ($number, $round_to, $split) = @_;

	my @numbers;
	for(my $i=1; $i<=$split; $i++) {
		my $number_split = int(($number / $split) * $i);
		my $rounded_number = &round_number_to_nearest_number($number_split, $round_to);
		$numbers[($i - 1)] = $rounded_number;
		#print "$i -1) $number (to nearest $round_to) = $rounded_number\n";
	}
	return \@numbers;
}

sub find_suitable_decanumber_to_round_to {
	my ($number, $split) = @_;
	my $smallest_part = int($number / $split);

	# Round smallest part to decanumber 
	my $decanumber = (substr $smallest_part, 0, 1);
	my $zeros;
	for(my $i=1; $i<length($smallest_part); $i++) { $zeros .= "0"; }
	$decanumber .= $zeros;

	# Get values then check for redundancy
	my $values = &split_and_round_numbers_to_nearest_number($number, $decanumber, $split);
	my %redundant;
	foreach(@{$values}) {
		if(defined $redundant{$_}) { warn "redundancy for $_ - must make new values!\n"; }
		$redundant{$_} = 1;
	}

	# This might be useful if I need to remove redundant numbers (i.e. inc/decrement by decanumber)
	#my $decanumber_increment = ('1' . $zeros);
	#warn "Decanumber = $decanumber - and increment by $decanumber_increment\n";

	return $decanumber;
}

sub split_and_round_to_suitable_number_and_display_suitable_unit {
	my ($number, $split) = @_;
	my $round = &find_suitable_decanumber_to_round_to($number, $split);
	my $values = &split_and_round_numbers_to_nearest_number($number, $round, $split);
	my ($new_values, $units) = &display_smaller_units_in_terms_of_nt($values, $round);
	return ($values, $new_values, $units);
}

sub display_smaller_units_in_terms_of_nt {
	my ($values_array, $rounded) = @_;

	my @values_in_smaller_units;
	my $units = "nt";

	# How many zeros
	my $how_many_zeros = ($rounded =~ tr/0//);
	#warn "how many zeros: $how_many_zeros\n";

	if($how_many_zeros >= 9) { $units = "Gb"; }
	if(($how_many_zeros < 9) && ($how_many_zeros >= 6)) { $units = "Mb"; }
	if(($how_many_zeros < 6) && ($how_many_zeros >= 3)) { $units = "Kb"; }

	#warn "i.e. $units\n";

	# Short hand
	if($units eq 'Gb') {
		foreach(@{$values_array}) {
			my $value = $_;
			$value =~ s/000000000//;
			push @values_in_smaller_units, $value;
		}
	}
	if($units eq 'Mb') {
		foreach(@{$values_array}) {
			my $value = $_;
			$value =~ s/000000//;
			push @values_in_smaller_units, $value;
		}
	}
	if($units eq 'Kb') {
		foreach(@{$values_array}) {
			my $value = $_;
			$value =~ s/000//;
			push @values_in_smaller_units, $value;
		}
	}

	return (\@values_in_smaller_units, $units);
}

1;
