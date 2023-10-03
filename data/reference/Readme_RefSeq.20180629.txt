6/29/2018 From Karsten Krug and Karl Clauser

This document describes the procedure to generate reference protein sequence databases for CPTAC 3.0

I. Content of this directory:
--------------------------
intermediate_2018-06-29_hg38					| intermediated files for hg38 assembly
---- 1_ucsc														| source files downloaded from UCSC table browser
---- 2_annotation_files_customProDB/	| annotation file created by customProDB
---- 3_annotate_fasta									| annotated protein sequence fasta file
---- 4_annotate_fasta_nr									| annotated protein sequence fasta file, with redundant sequences removed
intermediate_2018-06-29_mm10					| intermediated files for mm10 assembly, e.g. source files from UCSC table browser		
---- 1_ucsc														| source files downloaded from UCSC table browser
---- 2_annotation_files_customProDB/	| annotation file created by customProDB
---- 3_annotate_fasta									| annotated protein sequence fasta file
---- 4_annotate_fasta_nr									| annotated protein sequence fasta file, with redundant sequences removed
mitoContams							        			| FASTA files of lab contaminants and mitochondrial proteins
customProDBBI_1.13.3.tar.gz						| Broad version of customProDB R-package used to generate annotation files
Readme_RefSeq.20180629.txt						| this file


II. Final FASTA databases:
--------------------------
RefSeq.20180629_Human_Mouse_ucsc_hg38_mm10_cpdbnr_mito_264contams.fasta
RefSeq.20180629_Human_ucsc_hg38_cpdbnr_mito_264contams.fasta


III. Genomic annotation files: 
--------------------------
All proteins in the FASTA files (except contaminants and mitochondrial proteins) were mapped to the reference genome (hg38, mm10) and all 
annotation files were stored as RData-objects in corresponding subfolders, e.g  intermediate_2018-05-28_hg38/2_annotation_files_customProDB/.
These annotation files comprise id mapping tables (ids.RData), exon annotation (e.g. chr, strand, start and end positions) (exon_anno.RData), 
protein coding sequences (procodingseq.RData, proseq.RData) as well as a txdb.splite file containing all metadata in a single TxDb object 
which can be utilized by the GenomicFeatures R-package. The annotation files were created using the customProDB R-package, see below for a detailed description.


IV. Details about creation of these FASTA files is described below:
-----------------------------------------------------------
RefSeq.20180629_Human_ucsc_hg38_cpdbnr.fasta
RefSeq.20180629_Mouse_ucsc_mm10_cpdbnr.fasta

The description below is very thorough and it may help to start by answering the likely question, why didn't we just go to the the RefSeq website and download
their FASTA files?
The detail below was executed and documented with 3 major goals:
1) Maximize alignment to the genome for subsequent proteogenomic applications. Hence our FASTA files contain no RefSeq XP_ accession numbers.
2) Maximize annotation with standard HUGO gene symbols.
3) Keep a record so that 3 years from now we can explain how we obtained what we are using.

If some of the details like file names below seem not quite right it is probably a result of this document being updated from changes related to previous database generation efforts in
20160914 and 20171003.

FASTA files and accompanying annotation were downloaded from the UCSC table browser using hg19/hg38/mm10 genome assemblies. The annotation files were processed by the 'PrepareAnnotationRefseq'-function part of a mofified version of the CustomProDB R package (https://bioconductor.org/packages/release/bioc/html/customProDB.html).
The resulting 'ids.RData' file was then used to populate the FASTA headers with protein desciptions and gene symbols for ALL proteins in the
.fasta file. The overlap to HUGO was checked and is reported in venn diagrams. 

Detailed description (see R-vignette for customProDB package):
1) Downloading FASTA files form UCSC. 
The bullet below is taken from the R-vignette of CustomProDB and summarizes the steps to download coding sequence FASTA files:
- Go to UCSC Table Browser:	https://genome.ucsc.edu/cgi-bin/hgTables
- Choose genome:	Human or Mouse
- Choose assembly:	hg19 or hg38 or mm10
- Group:	Genes and Gene Prediction Tracks
- Track:	RefSeq Genes
- Table:	refGene
- Region:	genome
- Output format:	sequence
- Choose a filename:	e.g. 'ucsc_refseq_hg38_CDS_20180629.fasta' or 'ucsc_refseq_mm10_CDS_20180629.fasta'
- Press 'get output' button
- Then choose 'genomic' | CDS exons | one FASTA record per gene
- Press 'get sequence' button
Downloading the protein and mRNA sequence FASTA file is the same as above, just choose 'protein/mRNA' instead of 'genomic' after clicking the 'get output' button.

2) Preparing annotation files using CustomProDBBI v1.13.3
The two FASTA files (protein, CDS) downloaded in step 1) serve as input for CustomProDBBI:

>library(customProDBBI)
>PrepareAnnotationRefseq(genome='hg38', CDSfasta='ucsc_refseq_hg38_CDS_20180629.fasta', pepfasta='ucsc_refseq_hg38_protein_20180629.fasta', './annotation/', dbsnp=NULL, splice_matrix=T, COSMIC=F)

The function will generate several tables which get exported in folder './annotation' (.RData-files).

3) Extracting protein descriptions and gene symbols from annotation files generated in step 2) and generating the FASTA header:

>library(seqinr)
>fa <- read.fasta('ucsc_refseq_hg38_protein_20180629.fasta', seqtype='AA', as.string=F)
>load('./annotation/ids.RData')
>idx=which(ids$pro_name %in% sub('\\..*','',names(fa)))
>anno <- data.frame(matrix('', nrow=length(fa), ncol=ncol(ids), dimnames=list(sub("\\..*","",names(fa)), colnames(ids)) ))
>tmp <- data.frame(ids[idx, ])
>anno[match(tmp$pro_name, rownames(anno)), ] <- tmp
>anno <- data.frame(anno, ID=names(fa))
>fa.head <- apply(anno, 1, function(x) paste( x['ID'], ' ', ifelse(x['description'] != 'NA', x['description'], 'No annotation available'), ' ', 'GN=', ifelse(x['gene_name'] != 'NA', x['gene_name'], ''), ' [Homo sapiens]', sep='' ))
>write.fasta(fa, names=fa.head, file='ucsc_refseq_hg38_protein_20180629_customProDB.fasta', nbchar=50)


V. Appending lab contaminants and redundancy removal
------------------------------------------------------
To remove redundant FASTA entries the Spectrum Mill redundancy removal perl script was run on each of those databases. That script sorts all the proteins
in alphabetic order for accessions number, followed by descending order of sequence length and then removes entries that have 100% sequence identity
to an earlier entry in the FASTA file. The non-redundant human and mouse databases were then simply concatenated with each other, along with proteins
encoded by the mitochondrial genome, and common (non human/mouse)laboratory contaminants. See further detail below.

RefSeq.20180620_Human_mitochondria.fasta
RefSeq.20180620_Mouse_mitochondria.fasta
These files contain the 13 proteins for each species that are encoded in the mitochondrial genome and are often missing
when one obtains a reference sequence database based on assembly/mapping against the nuclear-encoded genome. These sequences
and accession numbers were obtained on 6/20/2018 from
	ftp://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/H_sapiens/protein/protein.fa.gz
	ftp://ftp.ncbi.nlm.nih.gov/refseq/M_musculus/mRNA_Prot/mouse.1.protein.faa.gz
The gene symbols were added manually and verified with
	human: https://www.genenames.org/cgi-bin/symbol_checker
	mouse: http://www.informatics.jax.org/batch

UniProtContams_264_20171227.fasta
A collection of 264 common labaoratory contaminants (non human/mouse) assembled from UniProt entries on 12/27/2017.
These are typically appended to all sequence databases at the Broad Insitute used for searching MS/MS spectra.
This includes a subset of 17 artifical entries (not in UniProt) accession numbers B99901-17 and  B99947-48 which represent
known artifacts in Promega trypsin and peptides routinely used for LC and MS calibration quality control.
