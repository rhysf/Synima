package rlines;
use strict;
use Exporter;
use Encode;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION = 0.1;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw();
%EXPORT_TAGS = (DEFAULT => [qw()], ALL =>[qw()]);
use FindBin qw($Bin);
use lib "$Bin/../modules";
use RoundNumber;
use Data::Dumper;

### rfarrer@broadinstitute.org

sub launch_R {
	my $R_file = $_[0];
	my $CMD = "Rscript $R_file";
	#my $CMD = "cat $R_file | R --vanilla";
	system($CMD);
	return 1;
}

sub close_device_print {
	my $fh = $_[0];
	print $fh "dev.off()\n";
	print $fh "q()\n\n";
	close $fh;
	return 1;
}

sub print_xlab_and_xaxis_for_windows {
	my ($xmax, $fh) = @_;
	my $number_of_points_on_axis = 9;
	my ($values, $rounded_values, $units) = roundnumber::split_and_round_to_suitable_number_and_display_suitable_unit($xmax, $number_of_points_on_axis);
	# units = "Mb for example", rounded_values = @(400, 800, 1200), values=@(4, 8, 12)

	# Create First point
	unshift @{$values}, '0';
	unshift @{$rounded_values}, '0';

	# Print Points of interest of xlab
	my $points_of_interest = "lablist.top<-as.vector(c(";
	foreach(@{$values}) { $points_of_interest .= "\"$_\", "; }
	$points_of_interest =~ s/,\s$//;
	#warn "points of interest = $points_of_interest))\n";
	print $fh "$points_of_interest))\n";

	# Axis
	print $fh "axis(1, at=lablist.top, line=0.4, lwd=2, labels=FALSE)\n";

	# Labels
	my $labels_of_interest = "axis(1, at=lablist.top, lab=c(";
	foreach(@{$rounded_values}) { $labels_of_interest .= "\"$_\", "; }
	$labels_of_interest =~ s/,\s$//;
	#warn "labels of interest = $labels_of_interest)\n";
	print $fh "$labels_of_interest), line=0.18, lwd=0, cex=1.5)\n";

	# Ylab
	print $fh "mtext(\"Position in genome ($units)\", outer = TRUE, side=1, line=2.4, cex=2)\n";
	return 1;
}

1;
