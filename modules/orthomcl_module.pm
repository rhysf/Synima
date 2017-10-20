#!/usr/bin/perl -w
## AUTHOR: Li Li, Feng Chen <fengchen@sas.upenn.edu>
## ORTHOMCL [2007-04-04] Version 1.4

## Copyright (C) 2004~2006 by University of Pennsylvania, Philadelphia, PA USA.
## All rights reserved.

## Before orthomcl.pl can be used, some variables (including directory variables
## or parameter variables) in orthomcl_module.pm need to be set, as described in
## ORTHOMCL_INSTALL.

package orthomcl_module;
use strict;
use Bio::SearchIO;
use Storable;
require Exporter;
use FindBin qw($Bin);
our $VERSION = 1.4;
our @ISA    = qw( Exporter );
our @EXPORT = qw(
                  executeFORMATDB
                  executeBLASTALL
                  executeMCL
                  constructDirectory
                  constructAllFasta
                  constructIDX_for_bpofile
                  constructSE_for_bpofile
                  retrieve_from_file
                  getline_from_bpofile
                  blast_parse
                  blastqueryab
                  read_ggfile
                  open_bpofile
                  nonredundant_list
                  numeric_value
                  write_matrix_index
                  mcl_backindex
                  write_parameter_log
                  write_endtime_in_parameter_log
                  dieWithUnexpectedError
                  write_log
                  write_bbh
                  mode5
                  $VERSION
                  $BLAST_PVALUE_CUTOFF_DEFAULT
                  $BLAST_FORMAT
                  $PERCENT_IDENTITY_CUTOFF_DEFAULT
                  $PERCENT_MATCH_CUTOFF_DEFAULT
                  $MAX_WEIGHT_DEFAULT
                  $MCL_INFLATION_DEFAULT
                  $all_fa_file
                  $genome_gene_file
                  $blast_file
                  $bpo_file
                  $bpo_idx_file
                  $bpo_se_file
                  $matrix_file
                  $index_file
                  $mcl_file
                  $mcl_bi_file
                  $parameter_log_file
                  $orthomcl_log_file
                  $bbh_file
                  @taxa
                  %gindex
                  %gindex2
                  %blastquery
                  %graph
                  %weight
);


# Software
our $BLASTALL                            = "blastall";
our $BLAST_FORMAT                        = "compact";  # "compact" corresponds to NCBI-BLAST's -m 8
                                                       # "full" corresponds to NCBI-BLAST's -m 0
                                                       # for WU-BLAST, make changes on subroutine executeBLASTALL
our $BLAST_NOCPU                         = 1;          # Useful when running BLAST on multi-processor machine
our $FORMATDB                            = "formatdb";

my $uname = $ENV{HOSTTYPE};
our $MCL                                 = "$Bin/mcl.$uname";

# path
our $PATH_TO_ORTHOMCL                    = "./";   # must end with "/" 
our $ORTHOMCL_DATA_DIR                   = $PATH_TO_ORTHOMCL . "/sample_data/";
our $ORTHOMCL_WORKING_DIR                = '';                            # will be changed with each run
our $ORTHOMCL_FORMER_RUN_DIR             = '';                            # if and only if user wants to use former run data
our $ORTHOMCL_TMP_DIR                    = '';                            # will be changed with each run

# cut-off variables
our $BLAST_PVALUE_CUTOFF_DEFAULT         = 1e-5;
our $PERCENT_IDENTITY_CUTOFF_DEFAULT     = 0;       # Both PercentIdentity and PercentMatch cutoff are 
our $PERCENT_MATCH_CUTOFF_DEFAULT        = 0;       # set to be 0 by default
our $MAX_WEIGHT_DEFAULT                  = 316;     # You need to see what's the second minimum p-value
                                                    # of your blast software. For example, if 9e-99 is 
                                                    # the case you should use -log(9e-99)=100.
our $MCL_INFLATION_DEFAULT               = 1.5;


# files
our $all_fa_file                         = "all.fa";
our $genome_gene_file                    = "all.gg";
our $blast_file                          = "all.blast";
our $bpo_file                            = "all.bpo";
our $bpo_idx_file                        = "all_bpo.idx";
our $bpo_se_file                         = "all_bpo.se";
our $matrix_file                         = "all_ortho.mtx";
our $index_file                          = "all_ortho.idx";
our $mcl_file                            = "all_ortho.mcl";
my  $bbh_file                            = "all_blast.bbh";
our $mcl_bi_file                         = "all_orthomcl.out";
our $parameter_log_file                  = "parameter.log";
our $orthomcl_log_file                   = "orthomcl.log";

# global variables storing data
our @taxa=();
our %gindex=();     # taxon_id -> reference to the array with gene_id as elements
our %gindex2=();    # gene_id -> taxon_id
our %blastquery=();

our %graph=();
our %weight=();

our $address_ref='';
our $blastquery_ref='';
our @mcl_index=();


# This subroutine is used to construct the directories to store
# the intermediate files and the final files.
# One or Two arguments:
# 1. String Variable: `date`
# 2. String Variable: (Optional) the existing directory, e.g. "July_21_2" where blast
#    file and bpo files are located. This option is useful with species
#    list unchanged (blast are already performed) but with orthomcl parameters
#    being changed.
# Last modified: 07/21/04
sub constructDirectory {
	my $date                 = $_[0];
	my $former_run_dir       = $_[1];

	$ORTHOMCL_WORKING_DIR = $PATH_TO_ORTHOMCL.(split(" ",$date))[1]."_".(split(" ",$date))[2];
	my $no=1;
	if (-e $ORTHOMCL_WORKING_DIR) {
		$no++;
		while (-e $ORTHOMCL_WORKING_DIR."_".$no) {
			$no++;
		}
		$ORTHOMCL_WORKING_DIR.="_$no";
	}
	$ORTHOMCL_WORKING_DIR.="/";
	system ("mkdir $ORTHOMCL_WORKING_DIR");
	$mcl_bi_file=$ORTHOMCL_WORKING_DIR.$mcl_bi_file;
	$parameter_log_file=$ORTHOMCL_WORKING_DIR.$parameter_log_file;
	$orthomcl_log_file=$ORTHOMCL_WORKING_DIR.$orthomcl_log_file;open (LOG,">$orthomcl_log_file") or die "can't create $orthomcl_log_file";
	$bbh_file=$ORTHOMCL_WORKING_DIR.$bbh_file;open (BBH,">$bbh_file") or die "can't create $bbh_file";
	write_log("### WORKING DIRECTORY: ###\n  $ORTHOMCL_WORKING_DIR\n\n");

	$ORTHOMCL_TMP_DIR=$ORTHOMCL_WORKING_DIR."tmp/";
	system ("mkdir $ORTHOMCL_TMP_DIR");
	
	if (defined $former_run_dir) {
		$ORTHOMCL_FORMER_RUN_DIR = $PATH_TO_ORTHOMCL.$former_run_dir."/tmp/";
		$all_fa_file             = $ORTHOMCL_FORMER_RUN_DIR.$all_fa_file;
		$genome_gene_file        = $ORTHOMCL_FORMER_RUN_DIR.$genome_gene_file;
		$blast_file              = $ORTHOMCL_FORMER_RUN_DIR.$blast_file;
		$bpo_file                = $ORTHOMCL_FORMER_RUN_DIR.$bpo_file;
		$bpo_idx_file            = $ORTHOMCL_FORMER_RUN_DIR.$bpo_idx_file;
		$bpo_se_file             = $ORTHOMCL_FORMER_RUN_DIR.$bpo_se_file;
	} else {
		$all_fa_file             = $ORTHOMCL_TMP_DIR.$all_fa_file;
		$genome_gene_file        = $ORTHOMCL_TMP_DIR.$genome_gene_file;
		$blast_file              = $ORTHOMCL_TMP_DIR.$blast_file;
		$bpo_file                = $ORTHOMCL_TMP_DIR.$bpo_file;
		$bpo_idx_file            = $ORTHOMCL_TMP_DIR.$bpo_idx_file;
		$bpo_se_file             = $ORTHOMCL_TMP_DIR.$bpo_se_file;
	}
	$matrix_file       = $ORTHOMCL_TMP_DIR.$matrix_file;
	$index_file        = $ORTHOMCL_TMP_DIR.$index_file;
	$mcl_file          = $ORTHOMCL_TMP_DIR.$mcl_file;
}

sub write_log {
	my $comment = $_[0];
	print LOG $comment;
	print STDERR $comment;
}
sub write_bbh {
	my $comment = $_[0];
	print BBH $comment;
}


# This subroutine is used to construct the directories to store
# the intermediate files and the final files.
# One or Two arguments:
# 1. String Variable: `date`
# 2. String Variable: the existing directory, e.g. "July_21_2" where matrix
#    file and index file are located. This option is useful with species
#    list changed but other orthomcl parameters unchanged.
# Last modified: 11/22/04
sub mode5 {
	my $date                 = $_[0];
	my $command              = $_[1];
	my $former_run_dir       = $_[2];
	my $usr_taxa_file        = $_[3];
	my $inflation            = $_[4];


	$ORTHOMCL_WORKING_DIR = $PATH_TO_ORTHOMCL.(split(" ",$date))[1]."_".(split(" ",$date))[2];
	my $no=1;
	if (-e $ORTHOMCL_WORKING_DIR) {
		$no++;
		while (-e $ORTHOMCL_WORKING_DIR."_".$no) {
			$no++;
		}
		$ORTHOMCL_WORKING_DIR.="_$no";
	}
	$ORTHOMCL_WORKING_DIR.="/";
	system ("mkdir $ORTHOMCL_WORKING_DIR");

	$ORTHOMCL_TMP_DIR=$ORTHOMCL_WORKING_DIR."tmp/";
	system ("mkdir $ORTHOMCL_TMP_DIR");
	
	my ($former_matrix_file,$former_index_file);
	if (defined $former_run_dir) {
		$ORTHOMCL_FORMER_RUN_DIR = $PATH_TO_ORTHOMCL.$former_run_dir.'/';

		open(FORMERLOG,$ORTHOMCL_FORMER_RUN_DIR.$parameter_log_file) || die;
		while (<FORMERLOG>) {
			if (/all\_fa\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$all_fa_file=$file;}
				else {$all_fa_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/genome\_gene\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$genome_gene_file=$file;}
				else {$genome_gene_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/blast\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$blast_file=$file;}
				else {$blast_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/blast\_parse\_out\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$bpo_file=$file;}
				else {$bpo_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/bpo\_idx\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$bpo_idx_file=$file;}
				else {$bpo_idx_file=$PATH_TO_ORTHOMCL.$file;}
			} elsif (/bpo\_se\_file\s+(\S+)/) {
				my $file=$1;
				if ($file=~/\//) {$bpo_se_file=$file;}
				else {$bpo_se_file=$PATH_TO_ORTHOMCL.$file;}
			}
		}
		$bbh_file                = $ORTHOMCL_FORMER_RUN_DIR.$bbh_file;
		$ORTHOMCL_FORMER_RUN_DIR.='tmp/';
		$former_matrix_file      = $ORTHOMCL_FORMER_RUN_DIR.$matrix_file;
		$former_index_file       = $ORTHOMCL_FORMER_RUN_DIR.$index_file;
	} else {die "wrong subroutine call!\n";}
	$matrix_file       = $ORTHOMCL_TMP_DIR.$matrix_file;
	$index_file        = $ORTHOMCL_TMP_DIR.$index_file;
	$mcl_file          = $ORTHOMCL_TMP_DIR.$mcl_file;

	$mcl_bi_file=$ORTHOMCL_WORKING_DIR.$mcl_bi_file;
	$parameter_log_file=$ORTHOMCL_WORKING_DIR.$parameter_log_file;
	$orthomcl_log_file=$ORTHOMCL_WORKING_DIR.$orthomcl_log_file;open (LOG,">$orthomcl_log_file") or die "can't open file";

	write_log("### WORKING DIRECTORY: ###\n  $ORTHOMCL_WORKING_DIR\n\n");

	read_ggfile($genome_gene_file);
	write_parameter_log($date,$command,5,"N/A","N/A","N/A",$inflation,"N/A");
	my %mcl_backindex;
	# read usr_taxa_file
	open (TAXA,$usr_taxa_file) or die "can't open $usr_taxa_file";
	while (<TAXA>) {
		if (/^([A-z]+)[\(\:\n]/) {
			my $taxon = $1;
			die "The taxon $1 is not defined!" unless (defined $gindex{$taxon});
			push(@mcl_index,@{$gindex{$taxon}});
		}
	}
	close(TAXA);
	for (my $i=0;$i<=$#mcl_index;$i++) {$mcl_backindex{$mcl_index[$i]} = $i;}

	# read former run index
	my @former_mcl_index;
	open (IDX,$former_index_file) or die "can't open $former_index_file";
	while (<IDX>) {
		if (/(\d+)	(\S+)/) {$former_mcl_index[$1]=$2;}
	}
	close(IDX);
	# read former run matrix
	my @matrix;
	open (MTX,$former_matrix_file) or die "can't open $former_matrix_file";
	while (<MTX>) {
		if (/^(\d+)\s+(.*)\$$/) {
			my $idx = $1;
			next unless (defined $mcl_backindex{$former_mcl_index[$idx]});
			my @score = split(" ",$2);
			foreach (@score) {
				if (/(\d+)\:(\S+)/) {
					next unless (defined $mcl_backindex{$former_mcl_index[$1]});
					push(@{$matrix[$mcl_backindex{$former_mcl_index[$idx]}]},$mcl_backindex{$former_mcl_index[$1]}.":".$2);
				}
			}
		}
	}
	close(MTX);

	my $size = scalar(@matrix);
	write_log("\nThere are $size sequences to cluster\n");
	open (MTX,">$matrix_file") or dieWithUnexpectedError("cannot write to file $matrix_file");
	print MTX "(mclheader\nmcltype matrix\ndimensions ".$size."x".$size."\n)\n\n(mclmatrix\nbegin\n\n";

	for (my $i = 0; $i <= $#matrix;$i++) {
		print MTX $i."    ";
		foreach (@{$matrix[$i]}) {
			print MTX $_." ";
		}
		print MTX "\$\n";
	}
	print MTX ")\n\n";
	close (MTX);

	write_log("Matrix($size X $size) file $matrix_file generated\n");

	open(IDX,">$index_file") or dieWithUnexpectedError("cannot write to file $index_file");
	for (my $i = 0;$i <= $#mcl_index; $i++) {
		print IDX "$i\t$mcl_index[$i]\n";
	}
	close(IDX);
	write_log("\nIndex file $index_file generated\n");
	
	executeMCL($matrix_file,$mcl_file,$inflation);
	mcl_backindex($mcl_file,$mcl_bi_file);
}


# Two arguments:
# 1. String Variable: Fasta file names, each representing one species, separated by comma, e.g. "Bme.fa,Ccr.fa,Eco.fa"
# 2. String Variable: Name of All_Fasta file
# Last modified: 10/17/07
sub constructAllFasta {
	my $fa_files=$_[0];
	my $all_fa_file=$_[1];
	my %seq_len;

	my @fafiles=split (",",$fa_files);
	my $totalgeneno=0;
	open (ALLFA,">$all_fa_file") or dieWithUnexpectedError("can't write to $all_fa_file!");
	foreach my $fafile (sort @fafiles) {
		open (FA,$ORTHOMCL_DATA_DIR.$fafile) or dieWithUnexpectedError("can't open $fafile!");
		my $taxon=(split (".fa",$fafile))[0];
		push (@taxa,$taxon);
		my $id;
		while (<FA>) {
			$_=~s/\r|\n//g;
			if (/^\>(\S+)/) {
				$id=$1;
				push (@{$gindex{$taxon}},$id);
				print ALLFA ">$id\n";
			} else {
				$seq_len{$id}+=length($_);
				print ALLFA "$_\n";
			}
		}
		close (FA);
		foreach (@{$gindex{$taxon}}) {$gindex2{$_}=$taxon;}
		$totalgeneno+=@{$gindex{$taxon}};
	}
	close (ALLFA);
	open (GG,">$genome_gene_file");
	foreach my $taxon (@taxa) {
		print GG "$taxon:";
		foreach (@{$gindex{$taxon}}) {
			print GG " $_";
		}
		print GG "\n";
	}
	close (GG);

	write_log("\nFASTA file <$all_fa_file> (".@fafiles." genomes, $totalgeneno sequences) generated!\n\n");
	return \%seq_len;
} ## constructAllFasta





# Two arguments:
# 1. String Variable: genome-gene file
# Last modified: 07/20/04
sub read_ggfile {
	my $gg_file=$_[0];
	my $totalgeneno=0;
	open (GG,$gg_file) or dieWithUnexpectedError("can't open $gg_file!");
	while (<GG>) {
		$_=~s/\r|\n//g;
		if (/(\S+)\(\d+\):(.*)/) {
			my $taxon=$1;
			my @genes=split (" ",$2);
			push (@taxa,$taxon);
			push (@{$gindex{$taxon}},@genes);
			foreach (@genes) {$gindex2{$_}=$taxon;}
			$totalgeneno+=scalar(@{$gindex{$taxon}});
		}
		elsif (/(\S+):(.*)/) {
			my $taxon=$1;
			my @genes=split (" ",$2);
			push (@taxa,$taxon);
			push (@{$gindex{$taxon}},@genes);
			foreach (@genes) {$gindex2{$_}=$taxon;}
			$totalgeneno+=scalar(@{$gindex{$taxon}});
		}
	}
	write_log("\nThere are ".@taxa." genomes, $totalgeneno sequences in $gg_file !\n\n");
	close (GG);
} ## read_ggfile




# One Argument:
# 1. String Variable: fasta file name
# -o t option is removed to solve gene id discrepancy
# Last modified: 05/10/06
sub executeFORMATDB {
	my $in = $_[0];
	write_log("\nRunning FormatDB\n");
	write_log("  Native FormatDB command:\n  $FORMATDB -i $in -p t\n\n");
	system ("$FORMATDB -i $in -p t");
} ## executeFORMATDB




# Four or Five Arguments:
# 1. String Variable: fasta file name
# 2. String Variable: blast out file name
# 3. String Variable: database name
# 4. String Variable: pvalue cutoff
# 5. (Optional) String Variable: "+number", parallel computing
# Last modified: 04/04/07
sub executeBLASTALL {
	my $in         = $_[0];
	my $out        = $_[1];
	my $db         = $_[2];
	my $pv_cutoff  = $_[3];
	my %blast_flag = %{$_[4]};

	write_log("\nRunning all-against-all BLAST\n");
	write_log("  Native BLASTALL command:\n  $BLASTALL -p blastp -i $in -d $db -e $pv_cutoff -o $out -m $blast_flag{'m'} -a $BLAST_NOCPU -v 1000 -b 1000\n");
	system ("$BLASTALL -p blastp -i $in -d $db -e $pv_cutoff -o $out -m $blast_flag{'m'} -a $BLAST_NOCPU -v 1000 -b 1000");
#-b and -v options are set to be the same large number
#-b Number of database sequences to show alignments for
#-v Number of database sequences to show one-line descriptions for
# by default, NCBI-BLAST sets v=500 and b=250; the difference in these
# two numbers can cause alignment information for the extra 250 hits
# unavailable (in NCBI-BLAST -m 0) for OrthoMCL, thus causing problems.
} ## executeBLASTALL




# Three Arguments:
# 1. String Variable: matrix file name
# 2. String Variable: mcl file name
# 3. Number Variable: inflation parameter for MCL algorithm
# Last modified: 07/19/04
sub executeMCL {
	my $matrix_file = $_[0];
	my $mcl_file    = $_[1];
	my $inflation  = $_[2];
	write_log("\nRun MCL program\n  $MCL $matrix_file -I $inflation -o $mcl_file\n");
	system("$MCL $matrix_file -I $inflation -o $mcl_file");
	if (-e $mcl_file) {
		write_log("\nMCL result $mcl_file generated!\n");}
	else {dieWithUnexpectedError("$mcl_file failed to be generated!");}
} ## executeMCL



# Three Arguments:
# 1. String Variable: blast out file
# 2. String Variable: blast parse out file
# 3. pv_cutoff
# Last modified: 05/10/06
sub blast_parse {
	my $blastfile        = $_[0];
	my $parseoutfile     = $_[1];
	my $pv_cutoff        = $_[2];
	my %seq_len          = %{$_[3]};
	my %blast_flag       = %{$_[4]};

	open (PARSEOUT,">$parseoutfile");
	write_log("\nParsing blast result!\n");
	my $searchio = Bio::SearchIO->new(-file   => $blastfile,
									  -format => $blast_flag{'format'}) or dieWithUnexpectedError("Blast parsing failed!");
	my $similarityid=1;
	while (my $result = $searchio->next_result()) {
		my $queryid=$result->query_name;
		my $querylen;
		if (defined $seq_len{$queryid}) {
			$querylen=$seq_len{$queryid};
		} else {
			$querylen=$result->query_length; # query length is not stored in BLAST m8 format, so querylen will be 0
		}
		while( my $hit = $result->next_hit ) {
			next unless numeric_pvalue($hit->significance) <= $pv_cutoff;
			my $subjectid=$hit->name;
			my $subjectlen;
			if (defined $seq_len{$subjectid}) {
				$subjectlen=$seq_len{$subjectid};
			} else {
				$subjectlen=$hit->length; # subject length is not stored in BLAST m8 format, so querylen will be 0
			}
			my $pvalue=numeric_pvalue($hit->significance);
			if ($blast_flag{'hsp'}) {
				my $simspanid=1;
				my $simspan='';
				my (@percentidentity,@hsplength);
				while( my $hsp = $hit->next_hsp ) {
					my $querystart=$hsp->start('query');
					my $queryend=$hsp->end('query');
					my $subjectstart=$hsp->start('sbjct');
					my $subjectend=$hsp->end('sbjct');
					$percentidentity[$simspanid]=$hsp->percent_identity;
					$hsplength[$simspanid]=$hsp->length('hit');
					$simspan.="$simspanid:$querystart-$queryend:$subjectstart-$subjectend.";
					$simspanid++;
				}
				my $sum_identical=0;
				my $sum_length=0;
				for (my $i=1;$i<$simspanid;$i++) {
					$sum_identical+=$percentidentity[$i]*$hsplength[$i];
					$sum_length+=$hsplength[$i];
				}
				my $percentIdent=int($sum_identical/$sum_length);
				print PARSEOUT "$similarityid;$queryid;$querylen;$subjectid;$subjectlen;$pvalue;$percentIdent;$simspan\n";
			} else {
				print PARSEOUT "$similarityid;$queryid;$querylen;$subjectid;$subjectlen;$pvalue;0;NULL\n";
			}
			$similarityid++;
		}
	}
	write_log("Parsing blast file finished\n");
	close(PARSEOUT);
} ## blast_parse




# This module is used to open BPO (Blast Parse Out) file, 
# to provide an filehandle for later use in blast query process.
# One Argument:
# 1. String Variable: BPO file name
# Last modified: 07/21/04
sub open_bpofile {
	my $bpo_file = $_[0];
	open (BPOFILE,$bpo_file) || dieWithUnexpectedError("can't open $bpo_file file"); # getline_from_bpofile subroutine 
																					 # will use this file handle
} # open_bpofile




# Make pvalue numeric, used by subroutine blast_parse
# One Arguments:
# 1. String Variable: pvalue
# Last modified: 07/19/04
sub numeric_pvalue {
	my $p=$_[0];
	if ($p=~/^e-(\d+)/) {return "1e-".$1;}
	else {return $p}
} # numeric_pvalue




# This subroutine is used to facilitate reading blast information
# from blastparseout file, instead of reading all blast result into
# memory.
# It works by storing every line (of blastparseout file)'s address
# into an array, and then uses Storable module to save this array
# into a file for future reference.
# Two Arguments:
# 1. String Variable: blast parse out file name
# 2. String Variable: index file name, storing index information of
#    blast parse out file.
# Last modified: 07/19/04
sub constructIDX_for_bpofile {
	my ($file,$idxfile)=@_;
	my @address;
	open (FILE,$file) or dieWithUnexpectedError("can't open $file file");
	push (@address, tell (FILE));
	while (<FILE>) {push (@address, tell (FILE));}
	store \@address, $idxfile;
	close (FILE);
} ## constructIDX_for_bpofile




# This subroutine is also used to facilitate reading blast information
# from blastparseout file, instead of reading all blast result into
# memory.
# It works by storing Starting and Ending line_ids of specific query_gene_id
# into a hash, and then uses Storable module to save this hash into a file
# for future reference.
# Two Arguments:
# 1. String Variable: blast parse out file name
# 2. String Variable: file name storing the hash
# Last modified: 07/19/04
sub constructSE_for_bpofile {
	my ($bpo_file,$bpo_se_file)=@_;
	my $lastquery='';my $lastsimid=0;
	open (BPOUT,$bpo_file) or dieWithUnexpectedError("can't open $bpo_file");
	while (<BPOUT>) {
		$_=~s/\r|\n//g;
		next unless ($_ ne '');
		my @bpout=split (";",$_);
		if ($lastquery ne $bpout[1]) {
			$blastquery{$bpout[1]}=$bpout[0];
			$blastquery{$lastquery}.=";".$lastsimid if $lastquery;
		}
		$lastquery=$bpout[1];
		$lastsimid=$bpout[0];
	}
	$blastquery{$lastquery}.=";".$lastsimid;
	store \%blastquery,$bpo_se_file;
	close(BPOUT);
} ## constructSE_for_bpofile




# This subroutine is used to retrieve @address and %blastquery
# from blast_idx_file and blast_SE_file
# Two arguments:
# 1. String Variable: blast_idx_file
# 2. String Variable: blast_SE_file
# Last modified: 07/19/04
sub retrieve_from_file {
	my $blastidxfile  = $_[0];
	my $blastsefile   = $_[1];
	$address_ref=retrieve ($blastidxfile);
	$blastquery_ref=retrieve ($blastsefile);
	%blastquery=%{$blastquery_ref};
} ## retrieve_from_file




# This subroutine is used to retrieve blast information from bpo file
# (blast parse out file), given the line_id.
# One Argument:
# 1. Reference Variable: reference to the address array
# 2. Number Variable: line_id of bpo file
# Last modified: 09/08/04
sub getline_from_bpofile {
	my $lineid         = $_[0];
	seek (BPOFILE,$address_ref->[$lineid-1],0);
	my $line=<BPOFILE>; 
	$line=~s/\r|\n//g;
	chop $line;
	my @bpout=split (";",$line);
	my ($pm,$pe);
	if ($bpout[5]==0) {$pm=0;$pe=0;}
#	elsif ($bpout[5]=~/(\d+)e(\-\d+)/) {$pm=$1;$pe=$2;}
	elsif ($bpout[5]=~/(\S+)e(\-\S+)/) {$pm=$1;$pe=$2;}  #For WU-BLAST their p-value has the pattern \d\.\de\-\d+  OR 0.
	else {$pm=$bpout[5];$pe=0;}
	return ($bpout[0],$bpout[1],$bpout[2],$bpout[3],,$bpout[4],$pm,$pe,$bpout[6],$bpout[7]);
} ## getline_from_bpofile




# This subroutine is used to retrieve blast information from bpo file
# (blast parse out file), given the query gene id and the subject gene
# id. If such information doesn't exist, zero is returned.
# Four Arguments:
# 1. String Variable: query gene id
# 2. String Variable: subject gene id
# Last modified: 07/19/04
sub blastqueryab {
	my $a              = $_[0];
	my $b              = $_[1];
	my ($s,$e);
	if (defined $blastquery_ref->{$a}) {
		($s,$e)=split (";",$blastquery_ref->{$a});
	} else {return 0;}
    foreach my $i ($s..$e) {
		my @c=(&getline_from_bpofile($i))[0,1,3,5,6,7];
		if ($c[2] eq $b) {return(@c);}
	}
	return 0;
} ## blastqueryab





# This subroutine is used to make a nonredundant list.
# One Argument:
# 1. Reference Variable: reference to an array
# Last modified: 07/19/04
sub nonredundant_list {
	my $list_ref=$_[0];
	my %nr=();
	foreach (@{$list_ref}) {$nr{$_}=1;}
	my @nr=sort (keys %nr);
	return \@nr;
} ## nonredundant_list





# This subroutine is used to generate input file (matrix file)
# for MCL and index file.
# Four arguments:
# 1. String Variable: matrix file name
# 2. String Variable: index file name
# Last modified: 02/24/05
sub write_matrix_index {

	my $matrix_file  = $_[0];
	my $index_file   = $_[1];

	my $size = scalar(keys %graph);
	write_log("\nThere are $size sequences to cluster\n");
	open (MTX,">$matrix_file") or dieWithUnexpectedError("cannot write to file $matrix_file");
	print MTX "(mclheader\nmcltype matrix\ndimensions ".$size."x".$size."\n)\n\n(mclmatrix\nbegin\n\n";

	my $i=0;
	my %mcl_index2;
	foreach my $p (sort keys %graph) {
		$mcl_index2{$p}=$i;$mcl_index[$i]=$p;
		$i++;
	}
	foreach my $p (sort keys %graph) {

		print MTX $mcl_index2{$p}."    ";
		foreach my $m (@{$graph{$p}}) {
			print MTX $mcl_index2{$m}.":".$weight{$p.' '.$m}." ";
		}
		print MTX "\$\n";
	}
	print MTX ")\n\n";

	close (MTX);

	write_log("Matrix($size X $size) file $matrix_file generated\n");

	#################################WRITE INDEX FILE######################################
	open(IDX,">$index_file") or dieWithUnexpectedError("cannot write to file $index_file");
	foreach my $id (sort { $mcl_index2{$a} <=> $mcl_index2{$b} } keys %mcl_index2) {
		print IDX "$mcl_index2{$id}\t$id\n";
	}
	close(IDX);
	write_log("\nIndex file $index_file generated\n");
}


# This subroutine is used to back index all gene_ids present 
# in MCL output and generate the final result.
# Two arguments:
# 1. String Variable: mcl file name
# 2. String Variable: mcl back_index file name
# Last modified: 11/13/06
sub mcl_backindex {
	
	my $mcl_file      = $_[0];
	my $mcl_bi_file   = $_[1];

	open (MCL,$mcl_file) or dieWithUnexpectedError("can't open $mcl_file");
	my $last=0;
	my $lastline='';
	my @mcl;
	while (<MCL>) {
#		chomp;chop;  #bug, reported by Robson Francisco de Souza.
                     #this may result in wrong gene index by chopping the last digit of line
                     #for new format of MCL output. Not a bug for older versions, e.g. mcl-02-063
                     #replaced to the following line: substitute the '$' with nothing, and 
                     #end reading data in the end of data block, preventing loading comment
                     #rows.
		chomp; s/\$$//; last if (/^\)/ && scalar(@mcl));

		if (/^(\d+)\s+(.*)\$/) {
			$mcl[$last]=$lastline;
			$last=$1;$lastline=$2;
		}
		elsif (/^(\d+)\s+(.*)/) {
			$mcl[$last]=$lastline;
			$last=$1;$lastline=$2;
		}
		elsif (/^\s+/) {$lastline.=$_;}
	}
	$mcl[$last]=$lastline;
	close (MCL);


	open (MCLBI,">$mcl_bi_file") or dieWithUnexpectedError("can't write to $mcl_bi_file");
	my $orthomcl_cluster_id=1;
	foreach my $mcl_cluster_id (0..$last) {
		$mcl[$mcl_cluster_id]=~s/\s+/ /g;
		my @g=split (" ",$mcl[$mcl_cluster_id]);
		next unless (scalar(@g)>=2);
		my @taxa=();
		foreach (@g) {
			my $taxon=$gindex2{$mcl_index[$_]};
			my $presence=0;
			foreach (@taxa) {
				if ($_ eq $taxon) {$presence=1;}
			}
			push (@taxa,$taxon) unless ($presence);
		}

		print MCLBI "ORTHOMCL".$mcl_cluster_id."(".scalar(@g)." genes,".scalar(@taxa)." taxa):	";
		foreach (@g) {
			print MCLBI " $mcl_index[$_]($gindex2{$mcl_index[$_]})";
		}
		print MCLBI "\n";
#No Species Cutoff
	}

	write_log("\n\nFinal ORTHOMCL Result: $mcl_bi_file generated!!!\n\n");

} ## mcl_backindex




# This subroutine is used to print out the parameters for future reference
# Eight Arguments:
# 1. String Variable: Start Time
# 2. String Variable: Command Line
# 3. Number Variable: OrthoMCL Mode number
# 4. Number Variable: P-value Cutoff
# 5. Number Variable: Percent Identity Cutoff
# 6. Number Variable: Percent Match Cutoff
# 7. Number Variable: MCL Inflation Parameter
# 8. Number variable: Maximum Weight
# Last modified: 07/21/04
sub write_parameter_log {
	my $starttime        = $_[0];
	my $command          = $_[1];
	my $mode             = $_[2];
	my $pv_cutoff        = $_[3];
	my $pi_cutoff        = $_[4];
	my $pmatch_cutoff    = $_[5];
	my $inflation        = $_[6];
	my $maximum_weight   = $_[7];
	open (PLOG,">$parameter_log_file") or die "can't open file";
	print PLOG "########################COMMAND######################\n";
	print PLOG "$command\n";
	print PLOG "######################PARAMETERS#####################\n";
	printf PLOG "%25s  %-10s\n","OrthoMCL Mode",$mode;
	printf PLOG "%25s  %-10s\n","P-value Cut-off",$pv_cutoff;
	printf PLOG "%25s  %-10s\n","Percent Identity Cut-off",$pi_cutoff;
	printf PLOG "%25s  %-10s\n","Percent Match Cut-off",$pmatch_cutoff;
	printf PLOG "%25s  %-10s\n","MCL Inflation",$inflation;
	printf PLOG "%25s  %-10s\n","Maximum Weight",$maximum_weight;
	printf PLOG "#######################SPECIES########################\n";
	foreach my $taxon (@taxa) {
		printf PLOG "%25s  %-20s\n",$taxon,scalar(@{$gindex{$taxon}})." genes";
	}
	printf PLOG "########################FILES########################\n";
	printf PLOG "%25s  %-50s\n","all_fa_file",$all_fa_file;
	printf PLOG "%25s  %-50s\n","genome_gene_file",$genome_gene_file;
	printf PLOG "%25s  %-50s\n","blast_file",$blast_file;
	printf PLOG "%25s  %-50s\n","blast_parse_out_file",$bpo_file;
	printf PLOG "%25s  %-50s\n","bpo_idx_file",$bpo_idx_file;
	printf PLOG "%25s  %-50s\n","bpo_se_file",$bpo_se_file;
	printf PLOG "%25s  %-50s\n","matrix_file",$matrix_file;
	printf PLOG "%25s  %-50s\n","index_file",$index_file;
	printf PLOG "%25s  %-50s\n","mcl_file",$mcl_file;
	printf PLOG "%25s  %-50s\n","mcl_bi_file",$mcl_bi_file;
	printf PLOG "%25s  %-50s\n\n","parameter_log_file",$parameter_log_file;
	print PLOG "######################START TIME#####################\n";
	print PLOG "$starttime\n";
	close (PLOG);
} ## write_parameter_log


# This subroutine is used to print out the end time in parameter log file
# One Argument:
# 1. String Variable: End Time
# Last modified: 09/03/04
sub write_endtime_in_parameter_log {
	my $endtime        = $_[0];
	open (PLOG,">>$parameter_log_file") or die "can't open file";
	print PLOG "########################END TIME#####################\n";
	print PLOG "$endtime\n";
	close (PLOG);
} ## write_endtime_in_parameter_log

# Last modified: 07/19/04
sub dieWithUnexpectedError {
	my $text = $_[ 0 ];
	die( "\n\n$0:\nUnexpected error (should not have happened):\n$text\n$!\n\n" );
} ## dieWithUnexpectedError


1;
