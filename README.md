
<img src="https://github.com/rhysf/Synima/blob/master/resources/logo.jpg?raw=true" width="400" height="400" />


## Introduction

Synima (Synteny Imager) is an orthology prediction pipeline and synteny viewer. The key features are:

* Orthologous genes are infered by either reciprocal best hits (RBH) from BLAST, OrthoMCL or Orthofinder. 
* Synteny is determined using DAGchainer and plotted using R
* All prerequisite programs are bundled with Synima, with wrapper scripts provided in /util/. 
* Several scripts can optionally launch in parallel  via LSF, GridEngine or UGER.

## Documentation

All documentation for Synima can be found at https://github.com/rhysf/Synima

## Support

For issues, questions, comments or feature requests, please check or post to the issues tab on github: https://github.com/rhysf/Synima/issues

## Prerequisites

* Perl
* Bio-Perl
* Python
* R
* Legacy-BLAST or BLAST+

## Getting started / examples

* Download software and run example Synima plot

    git clone git@github.com:rhysf/Synima.git  
    cd Synima/examples  
    perl ../SynIma.pl -a Repo_spec.txt.dagchainer.aligncoords -b Repo_spec.txt.dagchainer.aligncoords.spans

* Download software, run orthology pipeline and generate Synima plot

    git clone https://github.com/rhysf/Synima.git  
    cd Synima/examples  
    perl ../util/Create_full_repo_sequence_databases.pl -r ./Repo_spec.txt  
    perl ../util/Blast_grid_all_vs_all.pl -r ./Repo_spec.txt  
    perl ../util/Blast_all_vs_all_repo_to_OrthoMCL.pl -r ./Repo_spec.txt  
    ALTERNATIVELY 1: ../util/Blast_all_vs_all_repo_to_RBH.pl -r ./Repo_spec.txt  
    ALTERNATIVELY 2: ../util/Blast_all_vs_all_repo_to_Orthofinder.pl -r ./Repo_spec.txt  
    perl ../util/Orthologs_to_summary.pl -o all_orthomcl.out  
    perl ../util/DAGchainer_from_gene_clusters.pl -r ./Repo_spec.txt -c GENE_CLUSTERS_SUMMARIES.OMCL/GENE_CLUSTERS_SUMMARIES.clusters  
    perl ../SynIma.pl -a Repo_spec.txt.dagchainer.aligncoords -b Repo_spec.txt.dagchainer.aligncoords.spans``

## Description of the pipeline (Creating a sequence database)

* Synima visualises the output files from DAGChainer (aligncoords and aligncoords.spans files), which are tab delimited text files detailing the coordinates of sub-genomic regions of 
synteny between two or more genomes. 
* Having cloned a local copy of all the code using git clone, and navigated to the examples sub-folder, the first step is to create a 'repo sequence database'. 
* Create_full_repo_sequence_databases.pl reads in a Repository specification file (example Repo_spec file provided in examples) and outputs two fasta files 
(Repo_spec.txt.all.cds and Repo_spec.txt.all.pep) which are merged from 
each of the genome folders and used later.
* This first step is the most tricky - requiring that IDs in the GFF match the FASTA files. Warnings will alert the user to what ID's are being matched, and how many are matching. This step may need to be re-run until the correct settings or formatted files have been used.

The Input Repo_spec files take the format of:

``//
Genome CNB2
Annotation CNB2_FINAL_CALLGENES_1
//
Genome Cryp_gatt_IND107_V2
Annotation Cryp_gatt_IND107_V2_FINAL_CALLGENES_1
//``

* For each genome ID (E.g. CNB2) listed in the Repo_spec file, a corresponding 
sub-folder with the same name must be present in the same directory as the 
Repo_spec. 
* In each genome folder, a genome.fasta must be present, and named 
[genome-id].genome.fa E.g. CNB2/CNB2.genome.fa. 
* In each genome folder, annotation files should also be present, and named according to the details in the Data_repo. e.g., 3 additional files in each genome directory:

``[annotation-id].annotation.gff3
[annotation-id].annotation.cds
[annotation-id].annotation.pep``

* The GFF should have gene or mRNA features, with identifiers that are also in the 
.cds and .pep FASTA files, which contain the Coding sequence (cds) and peptide 
(pep) sequences respectively. 
* GFF ID's can be listed as such (or any other way that is specified explicitly for Create_full_repo_sequence_databases.pl): 

``cgbd CNB_WM276_v2 mRNA 450360 453023 . + . ID=012346;Parent=012345;Name=actin ``


* Gene ID's in the two FASTA files can be a single word i.e. >CBBG_0001 or they 
can have multiple fields e.g. 

``>01 gene_id=02 locus=03 name="ATPS" genome=Esch_coli analysisRun=Esch_coli_Augustus``

* If your files are downloaded from Ensembl, then utility scripts Ensembl_feature_table_to_gff.pl and Ensembl_reformat_fasta.pl may help parsing into the correct formats.

* Note: Gene and contig names should be alphanumerical (i.e., avoid symbols such as '=').

## Description of the pipeline (Predicting orthologous genes)

* With the sequence database made, the second step is to run all vs all BLAST hits using Blast_grid_all_vs_all.pl.

* Peptide or nucleotide alignments are possible, although peptide is generally recommended.

* BLAST searches can take a long time (especially with many genomes, or many predicted gene or proteins. Therefore, the option of distributing jobs 
to a cluster via LSF, SGE and UGE is provided (if available). 

* This step will create folders in each of the genome folders called RBH_blast_[PEP/CDS]. This 
step requires BLAST+ (makeblastdb and blastn/p) or BLAST legacy (formatdb and blastall) 
to be in $PATH.

* Next run either OrthoMCL, ORthofinder or reciprocal best hits (RBH) on the BLAST output 
using Blast_all_vs_all_repo_to_OrthoMCL.pl, Blast_all_vs_all_repo_to_Orthofinder.pl or Blast_all_vs_all_repo_to_RBH.pl,
respectively. 

* This will create an OMCL_outdir or RBH_outdir, containing all_orthomcl.out or PEP.RBH.OrthoClusters. 

* RBH will likely be less accurate than OrthoMCL or Orthofinder, but OrthoMCL at least has a limited number of genomes/genes that can be compared 
due to memory constraints.

* Next, summarise the OrthoMCL output (OMCL_outdir/all_orthomcl.out), 
or RBH output (RBH_outdir/PEP.RBH.OrthoClusters) or Orthofinder output (Orthofinder_outdir/
Orthogroups.csv) using Orthologs_to_summary.pl. 

* This step will create ortholog predictions in the output folders GENE_CLUSTERS_SUMMARIES.OMCL or 
GENE_CLUSTERS_SUMMARIES.RBH or GENE_CLUSTERS_SUMMARIES.Orthofinder respectively.

* The output of this step will also include some summary plots of the orthologs identified, and useful files for phylogenetics etc.

## Description of the pipeline (Visualising synteny)

* Run DAGChainer on the ortholog summary using DAGchainer_from_gene_clusters.pl.

* The final step is to run SynIma.pl on the aligncoords and aligncoords.spans output from DAGChainer. 

## Description of the pipeline (refining synteny plot)

* Once you have identified orthologs with the previous steps 1-5, you can re-run 
only this step with updated parameters to generate new figures. 

* If Synima finds the config.txt file (generated from the first time run, and in the same folder as
the figure, by default SynIma-output/config.txt), it will run using the parameters 
specified in this file (rather than use any updated parameters on the command 
line). 

* Config.txt includes a number of parameters that can change the appearance or layout
of the figure. It is recommended plotting both chromosome/contig synteny (c) and gene synteny (g)
separately, as either can give greater clarity depending on the input. 

* By default, synteny is shown as a partially transparent (alpha factor 0.5) azure4, although this can be changed to 
any other R color (E.g. http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf). Due to the 
color transparency, overlapping synteny will appear shaded. 

## Individual script details

* SynIma and wrapper script parameters are shown below, followed by their default settings given in [].

* Create_full_repo_sequence_databases.pl 

``Parameters: -r ./Repo_spec.txt []
Optional:   -f Feature wanted from GFF [mRNA]
	    -s Seperator in GFF description for gene names (\" ; etc) [;]
	    -d GFF description part number with the parent/gene info [0]
	    -m Remove additional comments in column [Parent=]
	    -v Verbose mode (report IDs missing) (y/n) [n]
Notes:      Will copy all transcripts and specified features from GFF into primary fasta files``

* Blast_grid_all_vs_all.pl

``Parameters: -r ./Repo_spec.txt []
Optional:   -t Type of alignment. Can be either PEP (peptide) or CDS (Coding sequence) [PEP]
            -c Number of best matches to capture between species [5]
            -s Number of top hits to capture in self-searches for paralogs [1000]  
            -e E-value cutoff [1e-20]
            -o Blast cmds outfile [blast.$type.cmds]
            -g Run commands on the grid (y/n) [n]
            -p Platform (UGER, LSF, GridEngine) [UGER]
            -q Queue name [short]
Notes:      Blast needs to be in PATH. 
            If BLAST+ (formatdb and blastn/p) is in PATH, that will be used. 
            Otherwise, BLAST legacy (formatdb and blastall) needs to be in PATH.``

* Blast_all_vs_all_repo_to_OrthoMCL.pl

``Parameters: -r ./Repo_spec.txt []
Optional:   -t Type of alignment. Can be either PEP (peptide) or CDS (Coding sequence) [PEP]
            -o Out directory [OMCL_outdir]``

* Blast_all_vs_all_repo_to_Orthofinder.pl

``Parameters: -r ./Repo_spec.txt []
Optional:   -t Type of alignment. Can be either PEP (peptide) or CDS (Coding sequence) [PEP]
            -o Out directory [Orthofinder_outdir]``

* Blast_all_vs_all_repo_to_RBH.pl

``Parameters: -r ./Repo_spec.txt
Optional:   -t Type of alignment. Can be either PEP (peptide) or CDS (Coding sequence) (PEP/CDS) [PEP]
	    -o Out directory [RBH_outdir]
	    -g Run commands on the grid (y/n) [n]
	    -p Platform (UGER, LSF, GridEngine) [UGER]
	    -q Queue name [short]``

* Orthologs_to_summary.pl 
``Parameters: -o Ortholog file (E.g. PEP.RBH.OrthoClusters, all_orthomcl.out, Orthogroups.csv) []
Optional:   -t Type of clustering (OMCL, RBH, Orthofinder) [OMCL]
            -d Outdir from Blast_all_vs_all_repo_to_OrthoMCL (if used) [OMCL_outdir]
	    	-r Repo Spec [./Repo_spec.txt]
            -p Repo Spec Peptide file [./Repo_spec.txt.all.PEP]``

* DAGchainer_from_gene_clusters.pl 

``Parameters: -r ./Repo_spec.txt []
            -c Ortholog cluster data (E.g. ORTHOMCLBLASTFILE.clusters) []
Optional:   -z File containing a list of genomes to restrict the analysis to []
            -i Minimum number of paired genes in a single dagchain [4]
            -o Cmds outdir [dagchainer_rundir]
            -l Cmds outfile [cluster_cmds]
            -g Run commands on the grid (y/n) [n]
            -p Platform (UGER, LSF, GridEngine) [UGER]
            -q Queue (hour, short, long) [short]
	    	-v Verbose (y/n) [n]
Notes:      GFF specifications (-f, -s, -d, -m) need to be the same as specified during 
            Blast_all_vs_all_repo_to_OrthoMCL.pl or Blast_all_vs_all_repo_to_RBH.pl";``

* Ensembl_feature_table_to_gff.pl
``Parameters: Feature table
Output: GFF``

* Ensembl_reformat_fasta.pl
``Parameters: FASTA (coding sequence/cds)
Output: FASTA``

* SynIma.pl

``Parameters: -c	./Config.txt [$cwd/SynIma-output/config.txt]
            -a	./Aligncoords []
            -b	./Aligncoords.spans []
Optional:   -e	Genome FASTA filename extension (e.g. ./SynIma/genome1/genome1.genome.fa etc.) [genome.fa]
            -t	Aligncoords.spans 2 []
	    	-u	Aligncoords.spans 3 []
	    	-k	Gene IDs 1 (1 per line) []
	    	-l	Gene IDs 2 (1 per line) []
	    	-o	Gene IDs 3 (1 per line) []
	    	-r	Run full program (y) or just create config (n) [y]
	    	-v	Verbose output (y/n) [n]
Plot Opts:  -i	Width of figure in pixels [1100]
	    	-j	Height of figure in pixels (num of genomes * 100)
            -g	Fill in chromosome/contig synteny (c) or gene synteny (g) [c]
	    	-z	Plot individual genes (y/n) [n]
            -x	Order of genomes from bottom to top seperated by comma
	    	-n	Genome labels from bottom to top seperated by comma
	    	-w	number of lines for left hand margin [12]
Notes:      Config.txt will be made automatically if not present, and read automatically if it is.
            Config.txt specifies order of genomes, chromosomes, colours, and other plot options. 
	    	Config.txt can be manually edited after creation.
	    	Default genome labels will be as they appear in aligncoords
	    	Order of genomes must have names as they appear in aligncoords
	    	Aligncoords.spans and Gene ID files will be highlighted according to the config``


* Ortholog_dist_per_genome_summary.pl 

``Parameters: -r [cluster_dist_per_genome.txt]
Optional: -c	clusters_and_unique.wAnnots
          -l	Seperate file with info to join []
	  	  -d	If opt_l, then from which column do i look for id? [0]
          -p	Printing options (n=none, c=look up names in opt_c and print those lines, l=get gene names from opt_c, then look up in opt_l and print those lines) [n]
	  	  -z	If opt_p = c, then restrict to 1 description - E.g. missing_in_Fo5176_454 []
Notes: Output files go to -c (wAnnots) directory. Specifically, this script divides orthogroups into categories of interest, including single-copy orthologs, which can
be used to construct phylogenetic trees, or genes that are not found in 1+ genome assemblies.``
