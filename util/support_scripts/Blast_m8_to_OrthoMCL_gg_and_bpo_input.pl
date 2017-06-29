#!/usr/bin/perl -w
use strict;

# Opening commands
my $usage = "Usage: perl $0 <m8_file> <orgSubstringLength> <outPrefix>\n";
die $usage unless(@ARGV eq 3);
my ($m8file, $substringLen, $out_prefix) = @ARGV;

main: {
	my %org_to_accs;
	my %seen;
	my $sim_counter = 0;
	open (my $bpo_fh, ">$out_prefix.bpo") or die $!;
	open (my $fh, $m8file) or die "Error, cannot open file $m8file";
	M8: while (my $line=<$fh>) {
		chomp $line;
		my @x = split (/\t/, $line);
		my ($accA, $accB, $Evalue, $per_id, $lendA, $rendA, $lendB, $rendB) = ($x[0], $x[1], $x[10], $x[2], $x[6], $x[7], $x[8], $x[9]);
		
		my $token = join ("$;", $accA, $accB);
		next M8 if ($seen{$token});
		$seen{$token} = 1;

		my $orgA = substr($accA, 0, $substringLen);
		my $orgB = substr($accB, 0, $substringLen);

		$org_to_accs{$orgA}{$accA} = 1;
		$org_to_accs{$orgB}{$accB} = 1;

		$per_id = int($per_id + 0.5);

		$sim_counter++;
		print $bpo_fh "$sim_counter;$accA;0;$accB;0;$Evalue;$per_id;1:$lendA-$rendA:$lendB-$rendB\n";
	}
	close $fh;
	close $bpo_fh;

	# write the gene list:
	open (my $gg_fh, ">$out_prefix.gg") or die $!;
	foreach my $org (keys %org_to_accs) {
		my @genes = keys %{$org_to_accs{$org}};
		print $gg_fh "$org: " . join (" ", @genes) . "\n";
	}
	close $gg_fh;
	print "Done.  See $$.gg and $$.bpo\n\n";
	exit(0);
}
	
__END__

FORMAT of all.gg or "usr_gg_file"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ath: At1g01190 At1g01280 At1g04160 ...
Hsa: Hs10834998 Hs10835119 Hs10835271 ...
Sce: YAL029c YAR009c YAR010c YHR023w ...

Each line stands for each genome. Each line starts with genome name, followed by a 
colon ":", and then followed by all the gene id's separated by space key " ".

FORMAT of all.bpo or "usr_bpo_file"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
1;At1g01190;535;At1g01190;535;0.0;97;1:1-535:1-535.
2;At1g01190;535;At1g01280;510;2e-56;29;1:69-499:28-474.
3;At1g01190;535;At1g11600;510;1e-45;27;1:59-531:21-509.

Each line represents each query-subject similarity relation. And all the info is
separated by ";", which are, in order, similarity id, query id, query length, 
subject id, subject length, BLAST E-value, percent identity, HSP info (each HSP
is in the format of HSP_id:query_start-query_end:subject_start-subject_end. 
different HSP info are seperated by "." )
IMPORTANT: 1. Similarity ID represents BPO file line id, so it should start 
              from 1 for the first line, and be consecutive for the whole file.
           2. BPO file is a parsing result from BLAST, so for each query gene
              id, its hits can't be scattered in the file, but should be listed 
              in ajacent lines.
           3. For BLAST m8 format (i.e. $BLAST_FORMAT="compact"), sequence length
              information is not stored. So when running OrthoMCL in mode 3,
              the corresponding columns (i.e. query length, and subject length)
              will be 0. Please do not use Percent match cutoff in this case.
