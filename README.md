# SynStack

These scripts can be used to evaluate the genomic neighborhood of a gene of interest in multiple genomes.
 
Prerequesites:
a) unix systems.

b) A folder with GFF3-files. The filename must start with a four letter species identifier (e.g. Drosophila melanogaster = DMEL) followed by an underscore. The GFF-files must be compressed ( DMEL_filename.gff3.gz)

c) A folder with protein files in fasta format. Here we also need a four letter species identifier in the beginning of the file name. The sequence title must start with the species id as well ( e.g. >DMEL_originalTitle).
The aminoacid sequence must follow the header in the next line and must be concatenated to only one line.

d) BLAST ( with makeblastdb and blastp)

e) MCL (https://micans.org/mcl/)

f) Bioperl with the following modules: 
Bio::Graphics, Bio::SeqFeature::Generic, Bio::Graphics::Glyph::segments, Bio::Graphics::Glyph::heat_map, Bio::Graphics::Glyph::generic

g) Gene tree of the genes we are interested in in NEWICK format. Titles in the gene tree must be identical to the titles provided to the pipeline.


The Pipeline is split into 4 seperated steps:
1) gather information about a specified genomic region around a gene of interest.

2) extract the proteinsequences from this genomic neighborhood and compare them to eachother. Identify clusters of similar sequences.

3) simplify the gff information so that all genes of interest have the same orientation, all genes have the same size.

4) visualise the extracted information so that all genes of interest are aligned in the and have the same color. each cluster of similar proteins have the same color and the cluster number is written on the gene. The order of the graphical output is given by the genetree to improve interpretation.

