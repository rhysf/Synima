#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $sequence_file = shift or die "Usage: $0 <sequence file> > out.sequence\n" ;
my $inseq = Bio::SeqIO->new('-file' => "<$sequence_file", '-format' => 'fasta' );
while (my $seq_obj = $inseq->next_seq ) {
    my $id = $seq_obj->id;
    my $seq = $seq_obj->seq;
    my $desc = $seq_obj->description;

    if ($desc =~ m/protein_id=([\w\d_\.]+)/) {
		my $protein_id = $1;
		#warn "$protein_id\n";
		print ">$protein_id $desc\n$seq\n";
    }
}

