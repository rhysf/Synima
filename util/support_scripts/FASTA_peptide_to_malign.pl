#!/usr/bin/env perl -w
use strict;
use Getopt::Std;
use FindBin qw($Bin);
use lib "$Bin/../../modules/";
use DataSpecFileParser;
#use Fasta_reader;
use read_FASTA;
use Run_Bsub;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case bundling);

# Opening statement
my $usage = "Usage: perl $0 -s <protein sequence to align file>\n
Optional: -c Use codon alignments instead of protein alignments
          -t Make tree (y/n) [y]\n";
our($opt_c, $opt_s, $opt_t);
getopt('cst');
die $usage unless ($opt_s);
if(!defined $opt_c) { $opt_c = 0; }
if(!defined $opt_t) { $opt_t = 'y'; }

# support scripts
my $make_codon_alignment_from_prot_alignment = "$FindBin::Bin/FASTA_make_codon_alignment_from_prot_alignment.pl";
#foreach($muscle_prog, $fastTree_prog, $make_codon_alignment_from_prot_alignment) { die "Cannot find $_ : $!\n" unless (-e $_); }
foreach($make_codon_alignment_from_prot_alignment) { die "Cannot find $_ : $!\n" unless (-e $_); }

main: {

	# Run MUSCLE alignment
	#my $cmd = "$muscle_prog -in $opt_s -out $opt_s.mfa";
	#my $cmd = "muscle -in $opt_s -out $opt_s.mfa";
	# version 5.1 syntax
	my $cmd = "muscle -align $opt_s -output $opt_s.mfa";
	if (-s "$opt_s.mfa") { warn "-malign file $opt_s.mfa already exists, not overwriting it\n"; }
        else { &process_cmd($cmd); }

	# Make Fast Tree
	my $alignment_file = "$opt_s.mfa";
        if($opt_t eq 'y') {
        	my $tree_file = "$opt_s.mfa.tree";
		my $cmd2 = "FastTree $opt_s.mfa > $tree_file";
        	if (-s $tree_file) { warn "-tree file $tree_file already exists, not over-writing it.\n"; }
        	else { &process_cmd($cmd2); }
	}

        # Convert to codon alignments
        if ($opt_c) {

		# Input/Output files
		my $cds_file = $alignment_file;
        $cds_file =~ s/pep\.mfa/cds/;
		my $codon_alignment_file = $alignment_file;
        $codon_alignment_file =~ s/\.pep\.mfa/\.codon\.mfa/;

		# Make new alignment from CDS
        my $cmd3 = "perl $make_codon_alignment_from_prot_alignment $alignment_file $cds_file > $codon_alignment_file";
		if (-s $codon_alignment_file) { warn "-not overwriting $codon_alignment_file\n"; }
		else { &process_cmd($cmd3); }

		# Buld a new tree
		if($opt_t eq 'y') {
			#my $cmd4 = "$fastTree_prog -nt -boot 0 $codon_alignment_file > $codon_alignment_file.tree";
			my $cmd4 = "FastTree -nt -boot 0 $codon_alignment_file > $codon_alignment_file.tree";
            		if (-s "$codon_alignment_file.tree") { warn "-not overwriting $codon_alignment_file.tree\n"; }
			else { &process_cmd($cmd4); }
		}
	}
	exit(0);
}

sub process_cmd {
	my ($cmd) = @_;
	warn "CMD: $cmd\n";
	my $ret = system($cmd);
	if ($ret) {
		#print STDERR "Error, cmd: $cmd died with ret: $ret";    
		die "Error, cmd: $cmd died with ret: $ret";
	}
	return;
}
