suppressMessages(library('ggplot2'))
suppressMessages(library('ggpubr'))
suppressMessages(library('ggdendro'))
suppressMessages(library(RColorBrewer))
suppressMessages(library(EnhancedVolcano))
suppressMessages(library(fgsea))
suppressMessages(library(scales))

# --------------------------------------------------
# fGSEA
# --------------------------------------------------
runGSEA <- function(
    df,
    rank.idx,
    group,
    feat,
    group_id,
    feature_id,
    gene_id,
    gmts,
    filt.pval=F,
    padj.thresh=NULL,
    minSize=3,
    nperm=10000,
    seed=NULL,
    ...
    ){
    # Run fGSEA
    # -------------------
    # df: input DataFrame
    # rank.idx: ranking for ranked fgsea
    # group: group to subset
    # feat: feature to subset
    # group_id: column ID for group
    # feature_id: feature ID for group
    # gene_id: gene Symbol ID
    # gmts: list of input .gmt files to use
    # filt.pval: whether or not to filter file for significant p-values
    # padj.thresh: threshold to filter adj. p-values if filt.pval = TRUE
    # minSize: fGSEA input
    # nperm: fGSEA input
    # seed: random seed
    if(!is.null(seed)){
      set.seed(seed)
    }

    df.filt <- df[(df[,group_id]==group) & (df[,feature_id]==feat),]
    df.filt <- df.filt[order(-df.filt$abs),]
    df.filt <- df.filt[ !duplicated(df.filt[,gene_id]),]

    R <- df.filt[,rank.idx]
    names(R) <- df.filt[,gene_id]
    e.df <- fgsea(gmts, R, nperm, minSize=minSize, ...)
    e.df$feature <- feat
    e.df$id <- group

    if(filt.pval){
      if (dim(e.df[e.df$padj<padj.thresh])[1] == 0){
        return(e.df[e.df$padj<padj.thresh])
        }
    }

    return(e.df)
}

runAllGSEA <- function(
    df,
    gmts,
    group_id='max_id',
    weight_id='max_norm',
    feature_id='feature',
    gene_id='geneSymbol',
    seed=NULL,
    minSize=3,
    nperm=10000,
    ...
    ){
    # Run All GSEA
    # -------------------
    # df: input DataFrame
    # gmts: list of input .gmt files to use
    # group_id: group to perform fgsea in
    # weight_id: rank for gsea
    # feature_id: feature type to perform fgsea for
    # gene_id: gene symbol ID
    # seed: random seed

    if(!is.null(seed)){
      set.seed(seed)
    }

    e_df <- list()

    # Convert columns to character
    df[,gene_id] <- as.character(df[,gene_id])
    #df$X <- as.character(df$X)
    df[,feature_id] <- as.character(df[,feature_id])

    # For each group / feature, perform fGSEA
    c=1
    for (group_ in unique(df[,group_id])){
        for (feature_ in unique(df[,feature_id])){
            e_df[[c]] <- runGSEA(df, weight_id, group_, feature_, group_id, feature_id, gene_id, gmts, minSize=minSize, nperm=nperm, ...)
            c = c + 1
        }
    }

    # Join results
    e_df <- do.call("rbind", e_df)
    e_df$leadingEdge <- as.character(e_df$leadingEdge)
    return(e_df)
}

# --------------------------------------------------
# Plotting
# --------------------------------------------------
plotGSEA <- function(e_df, pval.thresh=0.1, filter=NULL, palette='RdBu', h=13, w=15, s_color='black', remove.not.sig=TRUE){
    if(!is.null(filter)){
        e_df <- dplyr::filter(e_df, grepl(filter, pathway))
    }

    e_df$sig <- e_df$padj<pval.thresh
    e_df$logpval <- -log10(e_df$padj)

    if(remove.not.sig){
        e_df <- e_df[e_df$pathway %in% e_df[e_df$sig,]$pathway,]
    }

    ### Order axis by dendrogram
    # Load data
    X <- e_df[,c('pathway','feature','id','NES')]
    X$comb <- paste(X$feature,X$id)
    X <- reshape(X[,c('pathway','comb','NES')], timevar='comb', idvar='pathway', direction='wide',)
    rownames(X) <- X$pathway
    X$pathway <- NULL

    X[is.na(X)] <- 0

    # Build the dendrogram
    dend <- as.dendrogram(hclust(d = dist(x = X)))
    dendro.plot <- ggdendrogram(dend,rotate = TRUE)

    # Use dendrogram order to order colomn
    order <- order.dendrogram(dend) # dendrogram order
    e_df$pathway <- factor(x = e_df$pathway, levels = unique(e_df$pathway)[order], ordered = TRUE)

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        e_df,
        x="feature",
        y="pathway",
        fill = "NES",
        size="logpval",
        color=ifelse(e_df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(e_df$NES)) * c(-1, 1))+
        labs(x="", y="", fill="Enrichment", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
        facet_grid(. ~ id, scales = "free", space = "free")

    return(p)
}

plotPTMGSEA_v1 <- function(e_df, pval.thresh=0.1, filter=NULL, palette='RdBu', h=13, w=15, s_color='black'){
    if(!is.null(filter)){
        e_df <- dplyr::filter(e_df, grepl(filter, pathway))
    }

    e_df$sig <- e_df$fdr.pvalue<pval.thresh
    e_df$logpval <- -log10(e_df$fdr.pvalue)
    e_df <- e_df[e_df$pathway %in% e_df[e_df$sig,]$pathway,]

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        e_df,
        x="id",
        y="pathway",
        fill = "NES",
        size="logpval",
        color=ifelse(e_df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(e_df$NES)) * c(-1, 1))+
        labs(x="", y="", fill="Enrichment", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    return(p)
}

plotPTMGSEA <- function(e_df, pval.thresh=0.1, filter=NULL, palette='RdBu', h=13, w=15, s_color='black'){
    if(!is.null(filter)){
        e_df <- dplyr::filter(e_df, grepl(filter, pathway))
    }

    e_df$feature <- as.factor(e_df$feature)

    e_df$sig <- e_df$fdr.pvalue<pval.thresh
    e_df$logpval <- -log10(e_df$fdr.pvalue)
    e_df <- e_df[e_df$pathway %in% e_df[e_df$sig,]$pathway,]

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        e_df,
        x="feature",
        y="pathway",
        fill = "NES",
        size="logpval",
        color=ifelse(e_df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(e_df$NES)) * c(-1, 1))+
        labs(x="", y="", fill="Enrichment", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=45, hjust=1)) +
	facet_grid(. ~ id, scales="free", space="free")
    return(p)
}

plotCiberGSEA <- function(e_df, pval.thresh=0.1, filter=NULL, palette='RdBu', h=13, w=15, s_color='black'){
    if(!is.null(filter)){
        e_df <- dplyr::filter(e_df, grepl(filter, pathway))
    }

    e_df$id <- as.factor(e_df$id)

    e_df$sig <- e_df$padj<pval.thresh
    e_df$logpval <- -log10(e_df$padj)
    e_df <- e_df[e_df$pathway %in% e_df[e_df$sig,]$pathway,]

    ### Balloonplot
    options(repr.plot.width=w, repr.plot.height=h)

    p <- ggballoonplot(
        e_df,
        x="id",
        y="pathway",
        fill = "NES",
        size="logpval",
        color=ifelse(e_df$sig==T, s_color, "lightgrey")
        ) +
        scale_fill_distiller(palette=palette, limit = max(abs(e_df$NES)) * c(-1, 1))+
        labs(x="", y="", fill="Enrichment", size="-log10 Adj. P-val") + theme_linedraw() +
        theme(axis.text.x=element_text(angle=0))
    return(p)
}

plotVolcano <- function(de.df, w=10, h=12, xlim=NA, ylim=NA, gene.id='gene_name', genePTM=NULL, geneSelect=NULL, fcCut=0.05, pCut=0.1, lab.pval.thresh=NULL, ...){
    options(repr.plot.width=w, repr.plot.height=h)

    de.df[,gene.id] <- as.character(de.df[,gene.id])

    keyvals <- ifelse(
    de.df$adj.P.Val > pCut, 'grey',
      ifelse(abs(de.df$logFC) < fcCut, 'grey',
          ifelse(
              de.df$logFC < 0, 'royalblue',
              'red3'
          )
        )

    )

    names(keyvals)[keyvals == 'grey'] <- 'NS'
    names(keyvals)[keyvals == 'red3'] <- 'Up'
    names(keyvals)[keyvals == 'royalblue'] <- 'Down'

    xmin <- NA
    xmax <- NA
    if (!is.na(xlim)) {
        xmin <- xlim[1]
        xmax <- xlim[2]
    }

    if(!is.null(genePTM)){
        selectLab <- unique(de.df[de.df$geneSymbol %in% intersect(genePTM,de.df$geneSymbol),]$geneSite)
        if(!is.null(lab.pval.thresh)){
            selectLab <- de.df[de.df$geneSite %in% selectLab & de.df$adj.P.Val<lab.pval.thresh,]$geneSite
        }
        EnhancedVolcano(de.df,
            lab = de.df$geneSite,
            x = 'logFC',
            y = 'adj.P.Val',
            pCutoff=pCut,
            FCcutoff=fcCut,
            xlim=c(xmin,xmax),
            ylim=c(0,ylim),
            col=c('grey', 'grey', 'grey', 'red3'),
            ylab = bquote(~-Log[10]~ 'Adj. P-val'),
            colCustom = keyvals,
            selectLab=selectLab,
            ...
        )
    } else{
        EnhancedVolcano(de.df,
            lab = de.df[,gene.id],
            x = 'logFC',
            y = 'adj.P.Val',
            pCutoff=pCut,
            FCcutoff=fcCut,
            xlim=c(xmin,xmax),
            ylim=c(0,ylim),
            col=c('grey', 'grey', 'grey', 'red3'),
            ylab = bquote(~-Log[10]~ 'Adj. P-val'),
            colCustom = keyvals,
            ...
        )
    }
}

getVarSite <- function(x){
    return(as.character(strsplit(gsub("\\[|\\]", "", as.character(x)),', ')[[1]][1]))
}

loadClin <- function(FILE){
    clin.df = read.table(FILE, sep='\t', header=T)
    clin.df$log_adj_pval <- -log(clin.df$p_val_adj)
    clin.df$sig <- clin.df$p_val_adj<0.1
    clin.df$id <- as.factor(clin.df$id)
    clin.df$X <- factor(clin.df$X, levels=clin.df$X[order(-clin.df$p_val_adj)])
    clin.df$odds_r[clin.df$odds_r == Inf] <- max(clin.df$odds_r[is.finite(clin.df$odds_r)])
    clin.df <- clin.df[which(clin.df$p_val_adj <= 0.1),]
    return(clin.df)
}

plotClinicalFisher <- function(DIR, FILE, title="", outFile=NULL){
    # Load grading results
    clin.df <- loadClin(file.path(DIR,FILE))
    clin.df$src <- title
    # Plot
    options(repr.plot.width=12, repr.plot.height=3.5)
    mycols <- c("red", "steelblue")
    p <- ggballoonplot(
        clin.df,
	x="id",
	y="X",
	fill="odds_r",
	size="log_adj_pval",
	color=ifelse(clin.df$sig==T, "black", "grey"),
	rotate.x.text=TRUE,
	) +
	scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=1) +
	labs(x="", y="", size="-log10 Adj. P-val", fill="Odds Ratio") + theme_linedraw() +
	theme(axis.text.x=element_text(angle=45, hjust=1), legend.box="horizontal") +
	facet_grid(.~src, scales="free", space="free")
    if (!is.null(outFile)) {
        ggsave(outFile, p, dpi=320)
    }
    return(p)
}
