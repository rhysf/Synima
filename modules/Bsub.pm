#!/usr/local/bin/perl
package main;
our $SEE;
package Bsub;
use strict;
use warnings;
use Cwd;
use File::Basename;

my $WAITTIME = 15;
my $RETVAL_BIN_SIZE = 1000;
#my $RESORT_TO_POLLING_TIME = 60*60; # 1 hour (60 min * 60 seconds)
my $RESORT_TO_POLLING_TIME = 15*60; # 15 minutes
my %JOB_ID_TO_PREVTIME;
my $MAX_NODES = 500;

sub new {
	my $packagename = shift;
	my $params_href = shift;
	umask(0000);
	
	# must specify
	# log_dir: area to store return values #optional
	# cmd_list_aref: ordered list of commands to process
	# max_nodes : default 100
	# cmds_per_node: 10

	my $cmds_list_aref = $params_href->{cmds} or die "No commands specified";
	my $max_nodes = $params_href->{max_nodes} || $MAX_NODES;
	my $cmds_per_node = $params_href->{cmds_per_node} || 10;
	my $log_dir = $params_href->{log_dir} || cwd();
	my $queue = $params_href->{queue} || "broad";
	my $memory = $params_href->{memory};
	my $mount_test = $params_href->{mount_test};

	unless ($log_dir) {
		# create logging area
		$log_dir = "/local/scratch/bsubtmp";
		unless (-d $log_dir) {
		  mkdir ($log_dir);
		}
	}
	$log_dir .= "/bsub.J$$.$ENV{HOSTNAME}.$$." . time();
	unless (-d $log_dir) {
		mkdir $log_dir or die "Error, cannot mkdir $log_dir";
	}

	# write commands listing:
	open (my $fh, ">$log_dir/cmds_list.txt") or die $!;
	my $index = 0;
	foreach my $cmd (@$cmds_list_aref) {
		print $fh "index($index)\t$cmd\n";
		$index++;
	}
	close $fh;

	# finish logdir setup and object creation.
	my $cmds_dir = "$log_dir/cmds";
	my $retvals_dir = "$log_dir/retvals";
	my $monitor_dir = "$log_dir/monitor";
	foreach my $dir ($cmds_dir, $retvals_dir, $monitor_dir) {
		mkdir $dir or die "Error, cannot mkdir $dir";
	}

	my $num_cmds = scalar (@$cmds_list_aref);
	my $self = {
		num_cmds => $num_cmds,
		cmds_list => $cmds_list_aref,
		log_dir => $log_dir,
		cmds_dir => $cmds_dir,
		retvals_dir => $retvals_dir,
		monitor_dir => $monitor_dir,
		max_nodes => $max_nodes,
		cmds_per_node => $cmds_per_node,
		retvalues => 0, #later set to array ref with retvalues for each command.
		queue => $queue,
		nodes_in_progress => {},
		memory => $memory,
		mount_test => $mount_test,
		job_id_to_cmd_indices => {},  # job_id => [1,2,3,4,...]  so we know which cmds correspond to each grid job identifier.
		job_id_to_submission_time => {},
	};

	bless ($self, $packagename);
	return ($self);
}

sub bsub_jobs {
	my $self = shift;

	$self->_write_pid_file();
	my $max_nodes = $self->{max_nodes};
	my $num_cmds = $self->{num_cmds};
	my $num_cmds_launched = 0;
	my $num_nodes_used = 0;

	while ($num_cmds_launched < $num_cmds) {
		$num_cmds_launched = $self->_submit_job($num_cmds_launched);
		print STDERR "\r  CMDS: $num_cmds_launched / $num_cmds  [$num_nodes_used/$max_nodes nodes in use]	";
		$num_nodes_used = $self->_get_num_nodes_used();
		if ($num_nodes_used >= $max_nodes) {
			my $num_nodes_finished = $self->_wait_for_completions();
			$num_nodes_used -= $num_nodes_finished;
		}
	}

	print STDERR "\n* All cmds submitted to grid.  Now waiting for them to finish.\n";
	## wait for rest to finish
	while (my $num_nodes_finished = $self->_wait_for_completions()) { 
		$num_nodes_used -= $num_nodes_finished;
		print STDERR "\r  CMDS: $num_cmds_launched / $num_cmds  [$num_nodes_used/$max_nodes nodes in use]	";
	};
	
	print STDERR "\n* All nodes completed.  Now auditing job completion status values\n";
	
	$self->_get_exit_values();

	my $retvals_aref = $self->{retvalues};
	
	my $num_successes = 0;
	my $num_failures = 0;
	my $num_unknown = 0;

	foreach my $retval (@$retvals_aref) {
		if ($retval =~ /\d+/) {
		  if ($retval == 0) {
			 $num_successes++;
		  } else {
			 $num_failures++;
		  }
		} else {
		  $num_unknown++;
		}
	}

	$self->_write_result_summary($num_successes, $num_failures, $num_unknown);
	
	if ($num_successes == $num_cmds) {
		print "All commands completed successfully.\n";
	} else {
		print "Failures encountered:\n"
		  . "num_success: $num_successes\tnum_fail: $num_failures\tnum_unknown: $num_unknown\n";
	}
	print "Finished.\n\n";
}

sub _submit_job {
	my $self = shift;
	my $num_cmds_launched = shift;
	
	my $num_cmds = $self->{num_cmds};
	my $log_dir = $self->{log_dir};
	my $retvals_dir = $self->{retvals_dir};
	my $cmds_dir = $self->{cmds_dir};
	my $cmds_per_node = $self->{cmds_per_node};
	my $cmds_list_aref = $self->{cmds_list};
	my $monitor_dir = $self->{monitor_dir};

	my $orig_num_cmds_launched = $num_cmds_launched;
	
	
	my $shell_script = "$cmds_dir/J$$.S${num_cmds_launched}.sh";
	open (my $fh, ">$shell_script") or die $!;
	print $fh "#!/bin/sh\n\n";
	
	&_write_minimal_environment($fh);
	
	my $num_cmds_written = 0;

	my $monitor_started = "$monitor_dir/$num_cmds_launched.started";
	my $monitor_finished = "$monitor_dir/$num_cmds_launched.finished";

	my @cmd_indices_prepped;
	
	while ($num_cmds_launched < $num_cmds && $num_cmds_written < $cmds_per_node) {
		my $next_cmd_index = $num_cmds_launched; #always one less than the current index
		my $cmd_string = $cmds_list_aref->[ $next_cmd_index ];
		
		push (@cmd_indices_prepped, $next_cmd_index);
		
		my $retval_bin = int($next_cmd_index / $RETVAL_BIN_SIZE);
		
		my $retval_subdir = "$retvals_dir/$retval_bin";
		unless (-d $retval_subdir) {
			mkdir $retval_subdir or die "Error, cannot mkdir $retval_subdir";
		}

		print $fh "## Command index $next_cmd_index\n"
		  . "touch $monitor_started\n"
		  . "$cmd_string\n"
		  . 'echo $? >> ' . "$retval_subdir/entry_$next_cmd_index.ret\n\n";
		
		$num_cmds_launched++;
		$num_cmds_written++;
	}
	
	print $fh "\n" 
		. "rm -f $monitor_started\n"
		. "touch $monitor_finished\n"
		. "\n" 
		. "exit 0\n\n";
	
	
	close $fh;
	chmod (0775, $shell_script);
	
	print "Submitting: $shell_script to bsub\n" if $SEE;
	
	my $queue = $self->{queue};

	my $script_basename = basename($shell_script);
	my $cmd = "bsub -q $queue -e $shell_script.stderr -o $shell_script.stdout ";
	
	if (my $memory = $self->{memory}) {
		$cmd .= " -R \"rusage[mem=$memory]\" ";
	}
	if (my $mount_test = $self->{mount_test}) {
		$cmd .= " -E \"/broad/tools/NoArch/pkgs/local/checkmount $mount_test && [ -e $mount_test ]\" ";
	}
	if (my $group = $self->{group}) {
		$cmd .= " -G $group ";
	}
	
	$cmd .= " $shell_script 2>&1 ";
	

	my $job_id_text = `$cmd`;
	# print STDERR "\n$job_id_text\n";
	
	my $ret = $?;
	if ($ret) {
		print STDERR "BSUB failed to accept job: $cmd\n (ret $ret)\n";
		
		unlink $shell_script; # cleanup, try again later

		sleep(2*60); # sleep 2 minutes for now.  Give the system time to recuperate if a problem exists
		return ($orig_num_cmds_launched);
		
	}

	else {
		
		$shell_script = basename($shell_script);
		open (my $logdir_jobsfh, ">>$log_dir/job_ids.txt") or die "Error, cannot open file $log_dir/job_ids.txt";
		## get the job ID and log it:
		if ($job_id_text =~ /Job \<(\d+)\>/) {
		  my $job_id = $1;
		  print $logdir_jobsfh "$job_id\t$shell_script\n";
		  my $monitor_href = $self->{nodes_in_progress};
		  $monitor_href->{$monitor_finished} = $job_id;
		  
		  $self->{job_id_to_cmd_indices}->{$job_id} = \@cmd_indices_prepped;
		  $self->{job_id_to_submission_time}->{$job_id} = time();
		  
		}
		else {

		  die "Fatal error, couldn't extract Job ID from submission text: $job_id_text"; 
		  
		}
		close $logdir_jobsfh;
		
		# sleep($WAITTIME); # wait just a short while to give the system a few seconds to act on the submitted jobs.
		return ($num_cmds_launched);
	}
	
}

sub _get_num_nodes_used {
	my $self = shift;
	my $num_nodes_used = scalar (keys %{$self->{nodes_in_progress}});
	print "Num nodes currently in use: $num_nodes_used\n" if $SEE;
	return ($num_nodes_used);
}

sub _get_exit_values {
	my $self = shift;
	my $num_cmds = $self->{num_cmds};
	my @retValues;
	
	#print "Processing $retvals_dir\n";
	for (my $i = 0; $i < $num_cmds; $i++) {

		my $retval_file = $self->_get_ret_filename($i);
		
		#print "file: $retval_file\n";
		if (-s $retval_file) {
		  open (my $fh, $retval_file) or die $!;
		  my $retval_string = <$fh>;
		  $retval_string =~ s/\s//g;
		  $retValues[$i] = $retval_string;
		  close $fh;
		} else {
		  $retValues[$i] = "FILE_NOT_EXISTS";
		}
	}
	$self->{retvalues} = \@retValues;
}

sub get_failed_cmds {
	my $self = shift;
	my $retvalues_aref = $self->{retvalues};
	my $cmds_list_aref = $self->{cmds_list};

	my @failed_cmds;
	for (my $i = 0; $i <= $#$retvalues_aref; $i++) {
		my $retval = $retvalues_aref->[$i];
		if ($retval) {
		  push (@failed_cmds, 
				{ cmd => $cmds_list_aref->[$i],
				ret => $retval,
			 } );
		}
	}
	return (@failed_cmds);
}

sub _wait_for_completions {
	my $self = shift;
	
	print "sub _wait_for_completions()\n" if $SEE;
	
	my $nodes_in_progress_href = $self->{nodes_in_progress};
	
	my $seen_finished = 0;

	my @done;
	while (! $seen_finished) {
			 
		## check to see if there are any jobs remaining:
		if ($self->_get_num_nodes_used() == 0) {
		  ## no jobs in the queue
		  print "no nodes in use; exiting wait.\n" if $SEE;
		  return (0);
		}
		
		## check for finished jobs
		foreach my $monitor_file (keys %$nodes_in_progress_href) {
		  if (-e $monitor_file) {
			 push (@done, $monitor_file);
			 $seen_finished = 1;
		  }
		  else {
			 ## try polling the grid directly based on the job id
			 my $job_id = $nodes_in_progress_href->{$monitor_file};
			 
			 my $time_launched = $self->{job_id_to_submission_time}->{$job_id};
			 my $current_time = time();
			 
			 ## see if an hour has passed
			 if ($current_time - $time_launched >= $RESORT_TO_POLLING_TIME) {
				## poll the system directly:
				if (! $self->_job_running_or_pending_on_grid($job_id)) {
					
					push (@done, $monitor_file);
					$seen_finished = 1;
					
				}
				else {
					## reset submission time to delay next polling time
					$self->{job_id_to_submission_time}->{$job_id} = time();
				}
				
			 }
		  }
		  
		}
		if ($seen_finished) {
		  foreach my $monitor_file (@done) {
			 my $job_id = $nodes_in_progress_href->{$monitor_file};
			 print "job[$job_id]: $monitor_file is finished.\n" if $SEE;
			 delete $nodes_in_progress_href->{$monitor_file}; #remove from queue
			 delete $self->{job_id_to_cmd_indices}->{$job_id};
			 delete $self->{job_id_to_submission_time}->{$job_id};
		  }
		  return (scalar (@done)); #num jobs completed
		} 
		else {
		  ## wait a while and check again
		  print "waiting for jobs to finish.\n" if $SEE;
		  sleep($WAITTIME);
		}
	}
}


sub _write_pid_file {
	my $self = shift;
	my $log_dir = $self->{log_dir};
	open (my $fh, ">$log_dir/$ENV{HOSTNAME}.pid") or die $!;
	print $fh $$;
	close $fh;
}

sub _write_result_summary {
	my ($self, $num_successes, $num_failures, $num_unknown) = @_;
	my $status = ($num_failures == 0 && $num_unknown == 0) ? "success" : "failure"; 

	$self->{status} = $status;
	$self->{num_failures} = $num_failures;
	$self->{num_successes} = $num_successes;
	$self->{num_unknown} = $num_unknown;

	my $log_dir = $self->{log_dir};
	open (my $fh, ">$log_dir/bsub.finished.$status") or die $!;
	print $fh "num_successes: $num_successes\n"
		. "num_failures: $num_failures\n"
		. "num_unknown: $num_unknown\n";
	close $fh;
	
}

sub clean_logs {
	my $self = shift;
	my $log_dir = $self->{log_dir};
	my $cmd = "rm -rf $log_dir";
	system $cmd;
	return ($?);
}

sub _write_minimal_environment {
	my ($fh) = @_;

	print $fh <<_EOFENV_;

## add any special environment settings

echo HOST: \$HOSTNAME
echo HOST: \$HOSTNAME >&2

_EOFENV_

;
	return;
}

sub _get_ret_filename {
	my $self = shift;
	my ($cmd_index) = @_;

	my $retvals_dir = $self->{retvals_dir};

	my $retval_bin = int ($cmd_index / $RETVAL_BIN_SIZE);
	my $retval_file = $retvals_dir . "/$retval_bin/entry_$cmd_index.ret";
	
	return($retval_file);
}

sub _job_running_or_pending_on_grid {
	my $self = shift;
	my ($job_id) = @_;
	
	if (time() - $self->{job_id_to_submission_time}->{$job_id} < $RESORT_TO_POLLING_TIME) {
		return("TOO_SOON");
	}
	

	# print STDERR "Polling grid to check status of job: $job_id\n";
	
	my $response = `bjobs $job_id`;
	#print STDERR "Response:\n$response\n";

	foreach my $line (split(/\n/, $response)) {
		my @x = split(/\s+/, $line);

		if ($x[0] eq $job_id) {
		  my $state = $x[2];
		  if ($state =~ /^(DONE|EXIT|UNKWN|ZOMBI)$/) {
			 # better or worse, consider it done
			 return(0);
		  }
		  else {
			 $self->{job_id_to_submission_time}->{$job_id} = time();
			 return($state);
		  }
		}
	}
	
	print STDERR "-no record of job_id $job_id, setting as state unknown\n";
	return undef; # no status info
}

1; #EOM
