package IniReader;
use strict;
use warnings;
use Carp;

sub new {
	my $packagename = shift;
	my ($filename) = @_;

	# section -> att = value
	my $self = { section_to_att_val => {}, };

	my $current_section = "";
	open my $fh, '<', $filename or confess "Error, cannot open file $filename\n";
	while (<$fh>) {
		next if (/^[\;\#]/);
		next unless (/\w/);
		if (/\[([^\]]+)\]/) {
			$current_section = $1;
			$current_section = &_trim_flank_ws($current_section);
		} elsif (/^(.*)=(.*)$/) {
			my $att = $1;
			my $val = $2;

			$att = &_trim_flank_ws($att);
			$val = &_trim_flank_ws($val);
			$val =~ s/[\'\"]//g;
		
			$self->{section_to_att_val}->{$current_section}->{$att} = $val;
		}
	}
	close $fh;
	bless ($self, $packagename);
	return($self);
}

sub get_section_headings {
	my $self = shift;
	my @section_headings = keys %{$self->{section_to_att_val}};
	return(@section_headings);
}

sub get_section_attributes {
	my $self = shift;
	my $section = shift;
	my @attributes = keys %{$self->{section_to_att_val}->{$section}};
	return(@attributes);
}

sub get_value {
	my $self = shift;
	my ($section, $attribute) = @_;
	return ($self->{section_to_att_val}->{$section}->{$attribute});
}

sub _trim_flank_ws {
	my ($string) = @_;
	$string =~ s/^\s+|\s+$//g;
	return($string);
}

1; #EOM
