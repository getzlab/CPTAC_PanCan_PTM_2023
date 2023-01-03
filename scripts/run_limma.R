#!/usr/bin/Rscript
source("/home/yakiyama/CPTAC_PanCan_2021/scripts/limma.R")
suppressMessages(library(argparse))

parser <- ArgumentParser(description="Run differential expression for CPTAC data.")
parser$add_argument("--labels", help="Labels file.", required=TRUE, type="character")
parser$add_argument("--label_id", help="Labels id.", required=TRUE, type="character")
parser$add_argument("--feature_maps", help="File mapping features to genes.", required=TRUE, type="character")
parser$add_argument("--proteome", help="Protein input matrix.", default=NULL, required=FALSE, type="character")
parser$add_argument("--phosphoproteome", help="Phosphoproteome input matrix.", default=NULL, required=FALSE, type="character")
parser$add_argument("--phosphoproteome_res", help="Phosphoproteome residuals input matrix.", default=NULL, required=FALSE, type="character")
parser$add_argument("--acetylome", help="Acetylome input matrix.", default=NULL, required=FALSE, type="character")
parser$add_argument("--acetylome_res", help="Acetylome residuals input matrix.", default=NULL, required=FALSE, type="character")
parser$add_argument("--transcriptome", help="Transcriptome input matrix.", default=NULL, required=FALSE, type="character")
parser$add_argument("-o", "--output", help="Output directory.", required=TRUE, type="character")
parser$add_argument("--covar", help="Covariates file.", required=FALSE, default=NULL, type="character")
parser$add_argument("--covar_to_use", help="Covariates.", required=FALSE, default=c("cohort"), type="character", nargs="+")
parser$add_argument("--subset", help="Cluster IDs to subset.", default=NULL, required=FALSE, type="character", nargs="+")
parser$add_argument("--minObs", help="Minimum number of samples with measurement of a feature in each comparison group.", default=10, required=FALSE, type="integer")
args <- parser$parse_args()

dir.create(args$output)
cat(paste("   * Creating output directory:", args$output, "\n", sep=' '))

# -------------------------------------------------------------
# Load Reference Files
# -------------------------------------------------------------
# Load Labels
cat("   * Loading labels\n")
labels.df <- read.table(args$labels, sep='\t', header=TRUE, row.name=1)
cat(paste("     *", dim(labels.df)[1], "samples labeled\n", sep=' '))

# Subset labels
if (!is.null(args$subset)){
    labels.df <- labels.df[labels.df[,args$label_id] %in% args$subset,]
    cat(paste("     *", dim(labels.df)[1], "samples subsetted\n", sep=' '))
}

# Subset covariates
if (!is.null(args$covar)){
  cat("   * Loading covariates\n")
  covar.df <- read.table(args$covar, sep='\t', header=TRUE, row.names=1)
  covar.df <- covar.df[rownames(labels.df),]
  covar.df <- covar.df[,args$covar_to_use, drop=FALSE]
} else{
  covar.df <- NULL
}

# Protein Map
cat("   * Loading protein map\n")
id.map.df <- read.delim(args$feature_maps, sep='\t', header=TRUE, row.names=1)
id.map.df <- id.map.df[,c("geneSymbol"), drop=FALSE]
colnames(id.map.df) <- c("Description")

# -------------------------------------------------------------
# Run DE for each Assay
# -------------------------------------------------------------
cat("   * Running differential expression\n")
if (!is.null(args$proteome)){
  proteome.mat <- read.table(args$proteome, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
  cat(paste("     *", dim(proteome.mat)[1], "proteins loaded\n", sep=' '))
  cat("   * Running Proteome DiffExp\n")
  proteome.de <- RunDiffExp(proteome.mat, labels.df, args$label_id, id.map.df,  covar.df=covar.df, minObs=args$minObs)
  write.table(proteome.de, file.path(args$output, "proteome_de.tsv"), sep='\t', quote=FALSE)
}

if (!is.null(args$phosphoproteome)){
  phosphoproteome.mat <- read.table(args$phosphoproteome, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
  cat(paste("     *", dim(phosphoproteome.mat)[1], "phospho-sites loaded\n", sep=' '))
  cat("   * Running Phosphoproteome DiffExp\n")
  phosphoproteome.de <- RunDiffExp(phosphoproteome.mat, labels.df, args$label_id, id.map.df,  covar.df=covar.df, minObs=args$minObs)
  write.table(phosphoproteome.de, file.path(args$output, "phosphoproteome_de.tsv"), sep='\t', quote=FALSE)
}

if (!is.null(args$phosphoproteome_res)){
   phosphoproteomeRes.mat <- read.table(args$phosphoproteome_res, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
   cat(paste("     *", dim(phosphoproteomeRes.mat)[1], "phospho-sites loaded\n", sep=' '))
   cat("   * Running Phosphoproteome Residuals DiffExp\n")
   phosphoproteomeRes.de <- RunDiffExp(phosphoproteomeRes.mat, labels.df, args$label_id, id.map.df,  covar.df=covar.df, minObs=args$minObs)
   write.table(phosphoproteomeRes.de, file.path(args$output, "phosphoproteome_res_de.tsv"), sep='\t', quote=FALSE)
}

if (!is.null(args$transcriptome)){
  transcriptome.mat <- read.table(args$transcriptome, header=TRUE, sep='\t', row.names=1)[,rownames(labels.df)]
  cat(paste("     *", dim(transcriptome.mat)[1], "genes loaded\n", sep=' '))
  cat("   * Running Transcriptome DiffExp\n")
  transcriptome.de <- RunDiffExp(transcriptome.mat, labels.df, args$label_id, id.map.df, covar.df=covar.df, continuous=FALSE)
  write.table(transcriptome.de, file.path(args$output, "transcriptome_de.tsv"), sep='\t', quote=FALSE)
}

if (!is.null(args$acetylome)){
  acetylome.mat <- read.table(args$acetylome, header=TRUE, sep='\t', row.names=1)
  cat(paste("     *", dim(acetylome.mat)[1], "acetyl-sites loaded\n", sep=' '))
  a.samples <- intersect(rownames(labels.df),colnames(acetylome.mat))

  acetylome.mat <- acetylome.mat[,a.samples]
  a.labels.df <- labels.df[a.samples,]

  if (!is.null(args$covar)){
    a.covar.df <- covar.df[a.samples,]
  } else {
    a.covar.df <- NULL
  }

  cat("   * Running Acetylome DiffExp\n")
  acetylome.de <- RunDiffExp(acetylome.mat, a.labels.df, args$label_id, id.map.df, covar.df=a.covar.df, minObs=args$minObs)
  write.table(acetylome.de, file.path(args$output, "acetylome_de.tsv"), sep='\t', quote=FALSE)
}

if (!is.null(args$acetylome_res)){
  acetylomeRes.mat <- read.table(args$acetylome_res, header=TRUE, sep='\t', row.names=1)
  cat(paste("     *", dim(acetylomeRes.mat)[1], "acetyl-sites loaded\n", sep=' '))
  a.samples <- intersect(rownames(labels.df),colnames(acetylomeRes.mat))

  acetylomeRes.mat <- acetylomeRes.mat[,a.samples]
  a.labels.df <- labels.df[a.samples,]

  if (!is.null(args$covar)){
    a.covar.df <- covar.df[a.samples,]
  }

  cat("   * Running Acetylome DiffExp\n")
  acetylomeRes.de <- RunDiffExp(acetylomeRes.mat, a.labels.df, args$label_id, id.map.df, covar.df=a.covar.df, minObs=args$minObs)
  write.table(acetylomeRes.de, file.path(args$output, "acetylome_res_de.tsv"), sep='\t', quote=FALSE)
}

cat(paste("\n   * Outputs saved to", args$output, "\n", sep=' '))
