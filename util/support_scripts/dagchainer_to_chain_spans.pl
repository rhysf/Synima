#!/usr/bin/env perl
use strict;
use warnings;

my @chains;
my $curr_chain = undef;

while (<>) {
	chomp;
	my $line = $_;
	if (/^\#/) {
		
		my @x = split (/\s+/, $line);
		my $accA = $x[2];
		my $accB = $x[4];
		my $orientation = ($line =~ /\(reverse\)/) ? '-' : '+';
		$line =~ /num aligned pairs: (\d+)/ or die "Error, cannot find num aligned pairs";
		my $num_aligned_pairs = $1;

		
		
		$curr_chain = {
			accA => $accA,
			accB => $accB,
			orient => $orientation,
			num_gene_pairs => $num_aligned_pairs,
			coordsA => [],
			coordsB => [],
			
		};
		
		if ($line =~ /local inversion/) {
			$curr_chain->{local_inversion_flag} = 1;
		}
		

		push (@chains, $curr_chain);
		
	}
	else {
		my ($orgA, $contigA, $geneA, $lendA, $rendA, $posA, $matches,
			$orgB, $contigB, $geneB, $lendB, $rendB, $posB,
			$evalue) = split (/\t/);

		push (@{$curr_chain->{coordsA}}, $lendA, $rendA);
		push (@{$curr_chain->{coordsB}}, $lendB, $rendB);

	}
	
}

foreach my $chain (reverse sort {$a->{num_gene_pairs}<=>$b->{num_gene_pairs}} @chains) {
	
	if ($chain->{local_inversion_flag}) { next; }
	
	my $accA = $chain->{accA};
	my $accB = $chain->{accB};
	my $orient = $chain->{orient};
	my $num_gene_pairs = $chain->{num_gene_pairs};

	my @coordsA = @{$chain->{coordsA}};
	@coordsA = sort {$a<=>$b} @coordsA;
	my $lendA = shift @coordsA;
	my $rendA = pop @coordsA;

	my $lengthA = $rendA - $lendA + 1;
	

	my @coordsB = @{$chain->{coordsB}};
	@coordsB = sort {$a<=>$b} @coordsB;
	my $lendB = shift @coordsB;
	my $rendB = pop @coordsB;

	my $lengthB = $rendB - $lendB + 1;

	my $length_ratio = sprintf ("%.2f", $lengthA / $lengthB);

	print "$accA\t$lendA-$rendA\t$lengthA\t$accB\t$lendB-$rendB\t$lengthB\t$orient\t$num_gene_pairs\t$length_ratio\n";
	
}
exit(0);
