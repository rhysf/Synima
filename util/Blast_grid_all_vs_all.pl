#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../modules/";
use read_Blast;
use Synima;

### rfarrer@broadinstitute.org

# Opening commands
my $usage = "Usage: perl $0 -r <Repo_spec>
Optional: -t Type of alignment (PEP/CDS) [PEP]
          -c Number of best matches to capture between species [5]   # only single best hit
	  -s Number of top hits to capture in self-searches for paralogs [1000]  
	  -e E-value cutoff [1e-20]
	  -o Blast cmds outfile [blast.\$type.cmds]
	  -g Run commands on the grid (y/n) [n]
	  -p Platform (UGER, LSF, GridEngine) [UGER]
	  -q Queue name [short]\n
Note: BLAST legacy (formatdb and blastall) need to be in PATH (use BLAST)\n";
our($opt_r, $opt_t, $opt_c, $opt_s, $opt_e, $opt_o, $opt_g, $opt_p, $opt_q);
getopt('rtcseogpq');
die $usage unless($opt_r);
if(!defined $opt_t) { $opt_t = 'PEP'; }
if(!defined $opt_c) { $opt_c = 5; }
if(!defined $opt_s) { $opt_s = 1000; }
if(!defined $opt_e) { $opt_e = "1e-20"; }
if(!defined $opt_o) { $opt_o = "blast.$opt_t.cmds"; }
if(!defined $opt_g) { $opt_g = 'n'; }
if(!defined $opt_p) { $opt_p = 'UGER'; }
if(!defined $opt_q) { $opt_q = 'short'; }
die "Cannot open $opt_r : $!\n" unless(-e $opt_r);
die "-g is not n, N, y or Y: $opt_g\n" if($opt_g !~ m/^(n|N|y|Y)$/);
die "-t is not PEP or CDS: $opt_t\n" if($opt_t !~ m/^(PEP|CDS)$/);
die "-p is not UGER, LSF or GridEngine: $opt_p\n" if($opt_p !~ m/^(UGER|LSF|GridEngine)$/);

# Script to run commands
my $Run_Commands_python = "$Bin/support_scripts/Run_cmds_on_grid.py";
foreach($Run_Commands_python) { die "Cannot find $_ : $!\n" unless(-e $_); }

# Make commands for all-vs-all-search
my @blast_cmds = blastfile::make_all_vs_all_blast_search_cmds_from_repo($opt_r, $opt_t, $opt_c, $opt_s, $opt_e, $opt_o);

# Run commands on grid
if(($opt_g =~ m/y|Y/) && (-s $opt_o)) {
	my $cmd = "$Run_Commands_python --platform $opt_p --queue $opt_q --mem 2 --throttle_nodes 95 --cmds_per_node 10 $opt_o 2>&1 | tee $opt_o.log";
	synima::process_cmd($cmd);
} elsif($opt_g =~ m/n|N/) {
	open my $fh, '<', $opt_o or die "Cannot open $opt_o : $!\n";
	while(my $line=<$fh>) {
		chomp $line;
		if($line =~ m/\w+/) {
			warn "$line\n";
			system($line);
		}
	}
	close $fh;
} else { warn "Unsure of setting for running commands on grid or not: $opt_g\n"; }

# Clean up
`rm $opt_o` if($opt_o);
