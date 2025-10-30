#!/usr/local/bin/perl
package Run_Bsub;
use strict;
use warnings;
use Bsub;
use vars qw ($log_dir $keep_logs_on_failure);
use Cwd;
use Carp;

# by default, this is computed below. Set here specifically as needed.
our $CMDS_PER_NODE;  
my $NUM_NODES = 1000;

BEGIN {
	$log_dir = cwd;
	$keep_logs_on_failure = 1;
}

my $queue = "week";
my $memory = 4;
my $mount_test = undef;
my $max_nodes = undef;
my $group = undef;

sub set_queue {
	my ($q) = @_;
	$queue = $q;
	return;
}

# Should be in G
sub set_memory {
	my ($m) = @_;
	$memory = $m;
}

sub set_mount_test {
	my ($mt) = @_;
	$mount_test = $mt;
}

sub set_max_nodes {
	my ($mn) = @_;
	$max_nodes = $mn;
	return;
}

sub set_log_dir {
	my ($log_dir_setting) = @_;
	confess "Error, cannot find log directory $log_dir_setting" unless (-d $log_dir_setting);
	$log_dir = $log_dir_setting;
	return;
}

sub set_group {
	my ($group_setting) = @_;
	$group = $group_setting;
	return;
}

sub run {
	my @cmds = @_;
	my $num_cmds = scalar @cmds;

	my $targeted_cmds_per_node = ($CMDS_PER_NODE) ? $CMDS_PER_NODE : (int($num_cmds / $NUM_NODES) + 1);
	
	if ($targeted_cmds_per_node > 500) {
		$targeted_cmds_per_node = 500;
	}

	#if ($CMDS_PER_NODE && $CMDS_PER_NODE < $targeted_cmds_per_node) {
	#	$targeted_cmds_per_node = $CMDS_PER_NODE;
	#}

	my $bsubber = new Bsub({cmds=>\@cmds,
				log_dir => $log_dir,
				cmds_per_node => $targeted_cmds_per_node,
				queue => $queue,
				memory => $memory,
				mount_test => $mount_test,
				max_nodes => $max_nodes,
				group => $group,
				}
				);
	
	$bsubber->bsub_jobs();
	my $total_cmds = scalar (@cmds);

	if(my @failed_cmds = $bsubber->get_failed_cmds()) {
		my $num_failed_cmds = scalar (@failed_cmds);
		my $msg = "Sorry, $num_failed_cmds of $total_cmds failed = " . ($num_failed_cmds / $total_cmds * 100) . " % failure.\n";
		print $msg;
		print STDERR $msg;

		# 
		my @failed_cmd_lines;
		open (my $failed_fh, ">failed_cmds.$$") or die $!;
		foreach my $failed_cmd (@failed_cmds) {
			my $cmd = $failed_cmd->{cmd};
			my $ret = $failed_cmd->{ret};
			print $failed_fh "$cmd\nRET($ret)\n\n";
			push (@failed_cmd_lines, $cmd);
		}
		close $failed_fh;
		
		$bsubber->clean_logs() unless $keep_logs_on_failure;
		print "View file \'failed_cmds.$$\' for list of commands that failed.\n\n";

		# at least one job failed.
		return (@failed_cmd_lines); 
	} else {
		$bsubber->clean_logs();
		print "All $total_cmds completed successfully.\n\n";
		return(); # all good.
	}
}

package main;
use strict;
use warnings;
use File::Basename;

#### Test routine, run by executing module directly:
if (basename($0) eq 'Run_Bsub.pm') {
	my $ret = &Run_Bsub::run("ls");
	if ($ret) {
	   die "Error, test failed.\n";
	}
	else {
	   print "Test ran successfully.\n";
	}
	exit($ret);
}

1;
