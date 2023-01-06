#!/usr/bin/Rscript

suppressMessages(library(xCell))
source('/home/sanand/cibersort/CIBERSORT.R')
source("../../funcs/funcsR.R")
suppressMessages(library("ConsensusClusterPlus"))
suppressMessages(library(readr))
suppressMessages(library(ImmuneSubtypeClassifier))
suppressMessages(library(estimate))

# --------------------------------------------------
# Paths
# --------------------------------------------------
cat("    * Loading Data\n")
TRANSCRIPTOME_PATH <- "output/rna_tpm_gene_name.tsv.gz"
PROTEOME_PATH <- "../../data/processed/061721_genecentric/raw/proteome_X.tsv.gz"
OUT_DIR <- "ImmuneDeconv"
CIBERSORT_SIG <- "/home/sanand/cibersort/CIBERSORT_files/LM22.txt"

# Load
tpm.df <- read.delim(TRANSCRIPTOME_PATH, sep='\t', row.names=1, header=T)
prot.df <- read.delim(PROTEOME_PATH, sep='\t', row.names=1, header=T)[,colnames(tpm.df)]

cat(paste("    * Creating output directory: ", OUT_DIR, "\n"), sep=' ')
dir.create(OUT_DIR)

# --------------------------------------------------
# Run xCell
# --------------------------------------------------
xCell.out <- file.path(OUT_DIR, "xCell")
dir.create(xCell.out)

cat("    * Runing xCell (Transcriptome)\n")
rna_xcell <- xCellAnalysis(tpm.df)
write.table(rna_xcell, file.path(xCell.out, "rna_xcell.tsv"), sep='\t')

cat("    * Runing xCell (Proteome)\n")
prot_xcell <- xCellAnalysis(prot.df)
write.table(prot_xcell, file.path(xCell.out, "prot_xcell.tsv"), sep='\t')

# --------------------------------------------------
# Run CIBERSORT
# --------------------------------------------------
CIBER.out <- file.path(OUT_DIR, "cibersort")
dir.create(CIBER.out)

cat("    * Runing CIBERSORT\n")
set.seed(0, kind = "L'Ecuyer-CMRG" );
cibersort.results <- CIBERSORT(CIBERSORT_SIG, TRANSCRIPTOME_PATH)
write.table(cibersort.results, file.path(CIBER.out, "cibersort.tsv"), quote=FALSE, sep='\t')

# --------------------------------------------------
# Run ESTIMATE
# --------------------------------------------------
ESTIMATE.out <- file.path(OUT_DIR, "estimate")
dir.create(ESTIMATE.out)

fCommonGenes <- function(df, output.f, id = c("GeneSymbol", "EntrezID")){
    merged.df <- merge(common_genes, df, by.x = id, by.y = "row.names")
    rownames(merged.df) <- merged.df$GeneSymbol
    merged.df <- merged.df[, -1:-ncol(common_genes)]
    print(sprintf("Merged dataset includes %d genes (%d mismatched).",
        nrow(merged.df), nrow(common_genes) - nrow(merged.df)))
    outputGCT(merged.df, output.f)
}

cat("    * Runing ESTIMATE\n")
fCommonGenes(tpm.df, file.path(ESTIMATE.out , "cptac_tpm.gct"), id="GeneSymbol")
estimateScore(file.path(ESTIMATE.out , "cptac_tpm.gct"), file.path(ESTIMATE.out , "cptac_tpm_score.gct"), platform="illumina")

# --------------------------------------------------
# Run ImmuneSubtypeClassifier
# --------------------------------------------------
ISC.out <- file.path(OUT_DIR, "immune_subtype_classifier")
dir.create(ISC.out)

cat("    * Runing ImmuneSubtypeClassifier\n")
X.input <- as.matrix(tpm.df)
ISC.calls <- callEnsemble(X=X.input, geneids='symbol')
write.table(ISC.calls, file.path(ISC.out,"isc_results.tsv"), sep='\t')
