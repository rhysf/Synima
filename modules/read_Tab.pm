package tabfile;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);

### rfarrer@broadinstitute.org

sub save_one_column {
	my ($file, $col1_n) = @_;
	my %interest;
	my $tally = 0;
	warn "Saving $file...\n";
	open IN, "<$file" or die "Cannot open $file: $!\n";
	LIST: while(my $line=<IN>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($col1) = ($bits[$col1_n]);
		die "Undefined columns $col1_n: $line\n" if(!defined $col1);

		# Save
		$interest{$col1} = 1;
		$tally++;
	}
	close IN;
	warn "$tally found in $file\n";
	warn scalar(keys(%interest)) . " unique whole entries saved\n";
	return(\%interest);
}

sub save_three_to_four_columns {
	my ($file, $col1_n, $col2_n, $col3_n, $col4_n) = @_;
	my %interest;
	my $tally = 0;
	open IN, "<$file" or die "Cannot open $file: $!\n";
	LIST: while(my $line=<IN>) {
		chomp $line;
		my @bits = split /\t/, $line;
		my ($col1, $col2, $col3, $col4) = ($bits[$col1_n], $bits[$col2_n], $bits[$col3_n], $bits[$col4_n]);
		die "Undefined columns $col1_n, $col2_n, $col3_n, $col4_n: $line\n" if((! $col1) || (! $col2) || (! $col3) || (! $col4));

		# Save
		$interest{$col1}{$col2}{$col3} = $col4;
		$tally++;
	}
	close IN;
	warn "$tally found in $file\n";
	return(\%interest);
}

1;
