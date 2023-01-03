#!/usr/bin/Rscript
source("../../limma.R")
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Run differential expression for CPTAC data.")
parser$add_argument("--labels", help="Labels file.", required=TRUE, type="character")
parser$add_argument("--label_id", help="Labels id.", required=TRUE, type="character")
parser$add_argument("--feature_maps", help="File mapping features to genes.", required=TRUE, type="character")
parser$add_argument("--proteome", help="Protein input matrix.", required=TRUE, type="character")
parser$add_argument("--phosphoproteome", help="Phosphoproteome input matrix.", required=TRUE, type="character")
parser$add_argument("--acetylome", help="Acetylome input matrix.", default=NULL, required=FALSE, type="character")
parser$add_argument("--transcriptome", help="Transcriptome input matrix.", required=TRUE, type="character")
parser$add_argument("-o", "--output", help="Output directory.", required=TRUE, type="character")
parser$add_argument("--subset", help="Cluster IDs to subset.", default=NULL, required=FALSE, type="character", nargs="+")
parser$add_argument("--sva", help="Use SVA for linear regression.", default=FALSE, required=FALSE, type="logical")
args <- parser$parse_args()

dir.create(args$output)
cat(paste("   * Creating output directory:", args$output, "\n", sep=' '))

# Load Labels
cat("   * Loading labels\n")
labels.df <- read.table(args$labels, sep='\t', header=TRUE, row.name=1)
cat(paste("     *", dim(labels.df)[1], "samples labeled\n", sep=' '))

# Subset labels
if (!is.null(args$subset)){
    labels.df <- labels.df[labels.df[,args$label_id] %in% args$subset,]
    cat(paste("     *", dim(labels.df)[1], "samples subsetted", sep=' '))
}

# Load Files
cat("   * Loading matrices\n")
proteome.mat <- read.table(args$proteome, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
cat(paste("     *", dim(proteome.mat)[1], "proteins loaded\n", sep=' '))

phosphoproteome.mat <- read.table(args$phosphoproteome, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
cat(paste("     *", dim(phosphoproteome.mat)[1], "phospho-sites loaded\n", sep=' '))

if (!is.null(args$acetylome)) {
   acetylome.mat <- read.table(args$acetylome, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
   cat(paste("     *", dim(acetylome.mat)[1], "acetyl-sites loaded\n", sep=' '))
}

transcriptome.mat <- read.table(args$transcriptome, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
cat(paste("     *", dim(transcriptome.mat)[1], "genes loaded\n", sep=' '))

# Protein Map
cat("   * Loading protein map\n")
id.map.df <- read.table(args$feature_maps, sep='\t', header=TRUE, row.names=1, fill = TRUE )
id.map.df <- id.map.df[,c("geneSymbol"), drop=FALSE]
colnames(id.map.df) <- c("Description")

# -----------------------------------
# Run Differential Expression Tests
# -----------------------------------
cat("   * Running Proteome DiffExp\n")
proteome.de <- RunDiffExp(proteome.mat, labels.df, args$label_id, id.map.df, args$sva)
write.table(proteome.de, file.path(args$output, "proteome_de.tsv"), sep='\t', quote=FALSE)

cat("   * Running Phosphoproteome DiffExp\n")
phosphoproteome.de <- RunDiffExp(phosphoproteome.mat, labels.df, args$label_id, id.map.df, args$sva)
write.table(phosphoproteome.de, file.path(args$output, "phosphoproteome_de.tsv"), sep='\t', quote=FALSE)

if (!is.null(args$acetylome)) {
   cat("   * Running Acetylome DiffExp\n")
   acetylome.de <- RunDiffExp(acetylome.mat, labels.df, args$label_id, id.map.df, args$sva)
   write.table(acetylome.de, file.path(args$output, "acetylome_de.tsv"), sep='\t', quote=FALSE)
}

cat("   * Running Transcriptome DiffExp\n")
transcriptome.de <- RunDiffExp(transcriptome.mat, labels.df, args$label_id, id.map.df, args$sva, continuous=FALSE)
write.table(transcriptome.de, file.path(args$output, "transcriptome_de.tsv"), sep='\t', quote=FALSE)

cat(paste("\n   * Outputs saved to", args$output, "\n", sep=' '))
