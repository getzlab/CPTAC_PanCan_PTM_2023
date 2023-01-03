suppressMessages(library("edgeR"))
suppressMessages(library("limma"))
suppressMessages(library("SmartSVA"))
suppressMessages(library("qvalue"))
suppressMessages(library("DESeq2"))

RunLimmaTest <- function(mat.df, var.df, covariates.df=NULL, genes.df=NULL, use.sva=TRUE, n.sv=NULL, return_mod=FALSE) {
    # Computes differential expression using a combination of SmartSVA and Limma T-Test
    # Modified for just purely continuous data
    #
    # Args:
    #   mat.df: data.frame with continuous signal (#expression x #samples)
    #   var: variable of interest
    #   genes.df: data.frame mapping gene IDs to gene names
    #   n.sv: number of surrogates for SVA. Automatically determined if NULL.
    #
    # Returns:
    #   limma t-test output augmented with Storey q-values
    #   SVA surrogate variables
    stopifnot(dim(unique(var.df))[1]>1 && dim(var.df)[2]==1)
    var <- colnames(var.df)[1]
    design <- cbind(1, var.df)

    # define model
    mod <- model.matrix(as.formula(paste0('~', var)), var.df)
    inGroup.mat <- mat.df[, which(mod[,2]==1)]
    outGroup.mat <- mat.df[, which(mod[,2]==0)]
    
    # split matrix to calculate proportion of missing per group

    if (use.sva) {
        cat("  * running SmartSVA\n")
        v <- as.matrix(mat.df)
        Y.r <- t(resid(lm(as.formula(paste0('t(v) ~ ', var)), data=var.df)))

        if (is.null(n.sv)){
          n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
        }

        cat(paste0("  * SVs: ", n.sv, "\n"))
        sv.obj <- smartsva.cpp(as.matrix(mat.df), mod, mod0=NULL, n.sv=n.sv, alpha=1, B=200, VERBOSE=TRUE)

        # update model to include SVs
        mod <- model.matrix(as.formula(paste0('~', var, '+sv.obj$sv')), var.df)  # intercept, var, SVs
    }

    if (!is.null(covariates.df)) {
        covar.var.df <- cbind(var.df, covariates.df)

        if (use.sva){
          mod <- model.matrix(as.formula(paste0('~', paste(paste(colnames(covar.var.df),collapse = '+'),"+sv.obj$sv"))), covar.var.df)
        } else{
          mod <- model.matrix(as.formula(paste0('~', paste(colnames(covar.var.df),collapse = '+'))), covar.var.df)
        }
    }

    cat(paste0("  * model matrix dimensions: ", dim(mod)[2], "\n"))

    if (return_mod) {
      return(mod)
    }

    # run limma
    fit <- lmFit(mat.df, design=mod, na.action=na.exclude)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=ncol(design), n=Inf, sort.by="none")

    # calculate q-values
    Q <- qvalue(res[, "P.Value"])
    res[, 'qval'] <- signif(Q$qvalue, 6)
    cat(paste0("Differentially expressed genes at 0.05 FDR: ", sum(res[, 'qval']<=0.05), "\n"))
    if (!is.null(genes.df)) {
        res[, 'gene_name'] <- genes.df[row.names(res), 'Description']
        res <- res[, colnames(res)[c(8,1:7)]]
    }

    #return(res)

    # Add missingness column to results file
    prop_missing <- rowSums(is.na(mat.df))/dim(mat.df)[2]
    res$propMissing <- prop_missing[rownames(res)]
    prop_missingIn <- rowSums(is.na(inGroup.mat))/dim(inGroup.mat)[2]
    prop_missingOut <- rowSums(is.na(outGroup.mat))/dim(outGroup.mat)[2]
    res$propMissingIn <- prop_missingIn[rownames(res)]
    res$propMissingOut <- prop_missingOut[rownames(res)]

    if (use.sva) {
        sva.df <- data.frame(sv.obj$sv)
        rownames(sva.df) <- colnames(mat.df)
        return(list("res"=res, "C"=sva.df))
    } else {
        return(res)
    }
}

RunDiffExprAnalysisLimma <- function(counts.df, var.df, covariates.df=NULL, genes.df=NULL, use.sva=TRUE, n.sv=NULL, return_mod=FALSE) {
    # Computes differential expression using a combination of SmartSVA and voom-limma
    #
    # Args:
    #   counts.df: data.frame with read counts (#genes x #samples)
    #   var: variable of interest
    #   genes.df: data.frame mapping gene IDs to gene names
    #   n.sv: number of surrogates for SVA. Automatically determined if NULL.
    #
    # Returns:
    #   voom-limma output augmented with Storey q-values
    #   SVA surrogate variables
    stopifnot(dim(unique(var.df))[1]>1 && dim(var.df)[2]==1)
    var <- colnames(var.df)[1]

    design <- cbind(1, var.df)

    # apply edgeR normalization (TMM) to counts
    dge <- DGEList(counts=counts.df)
    dge <- calcNormFactors(dge)

    # define model
    mod <- model.matrix(as.formula(paste0('~', var)), var.df)

    if (use.sva) {
        cat("  * running SmartSVA\n")
        v <- voom(dge, design=mod)  # run on transformed counts

        Y.r <- t(resid(lm(as.formula(paste0('t(v$E) ~ ', var)), data=var.df)))

        if (is.null(n.sv)){
          n.sv <- EstDimRMT(Y.r, FALSE)$dim + 1
        }

        cat(paste0("  * SVs: ", n.sv, "\n"))
        sv.obj <- smartsva.cpp(as.matrix(v$E), mod, mod0=NULL, n.sv=n.sv, alpha=1, B=200, VERBOSE=TRUE)

        # update model to include SVs
        mod <- model.matrix(as.formula(paste0('~', var, '+sv.obj$sv')), var.df)  # intercept, var, SVs
    }

    if (!is.null(covariates.df)) {
        covar.var.df <- cbind(var.df, covariates.df)

        if (use.sva){
          mod <- model.matrix(as.formula(paste0('~', paste(paste(colnames(covar.var.df),collapse = '+'),"+sv.obj$sv"))), covar.var.df)
        } else{
          mod <- model.matrix(as.formula(paste0('~', paste(colnames(covar.var.df),collapse = '+'))), covar.var.df)
        }
    }

    cat(paste0("  * model matrix dimensions: ", dim(mod)[2], "\n"))

    if (return_mod) {
      return(mod)
    }

    # run limma
    v <- voom(dge, design=mod)
    fit <- lmFit(v, design=mod, na.action=na.exclude)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=ncol(design), n=Inf, sort.by="none")

    # calculate q-values
    Q <- qvalue(res[, "P.Value"])
    res[, 'qval'] <- signif(Q$qvalue, 6)
    cat(paste0("Differentially expressed genes at 0.05 FDR: ", sum(res[, 'qval']<=0.05), "\n"))
    if (!is.null(genes.df)) {
        res[, 'gene_name'] <- genes.df[row.names(res), 'Description']
        res <- res[, colnames(res)[c(8,1:7)]]
    }

    if (use.sva) {
        sva.df <- data.frame(sv.obj$sv)
        rownames(sva.df) <- colnames(counts.df)
        return(list("res"=res, "C"=sva.df))
    } else {
        return(res)
    }
}

RunDiffExp <- function(mat.df, labels.df, id, map.df, covar.df=NULL, continuous=TRUE){
    # Generalizes DE Tests for each Cluster
    #
    # mat.df: matrix input (features x samples)
    # labels.df: labels dataframe with cluster IDs
    # id: label of column to perform DE on
    # map.df: mapping of features to genes
    # covar.df: covariates to include
    # continuous: whether to run Limma expecting counts (Voom-transform) or as continuous variables
    cluster_ids <- unique(labels.df[,id])
    deClusterResults <- vector(mode = "list", length=length(cluster_ids))

    for(i in 1:length(cluster_ids)){
        cat(paste(" * Cluster:", cluster_ids[i], sep=' '))

        labels.df$de <- labels.df[,id]==cluster_ids[i]

        # Run Limma on Continuous
        if(continuous){
            deClusterResults[[i]] <- RunLimmaTest(
                mat.df,
                labels.df[,c('de'), drop=FALSE],
                genes.df=map.df,
                use.sva=FALSE,
                covariates.df=covar.df
            )
        } else{
            deClusterResults[[i]] <- RunDiffExprAnalysisLimma(
                mat.df,
                labels.df[,c('de'), drop=FALSE],
                genes.df=map.df,
                use.sva=FALSE,
                covariates.df=covar.df
            )
        }

        deClusterResults[[i]]$id <- cluster_ids[i]
        deClusterResults[[i]]$index <- rownames(deClusterResults[[i]])
        rownames(deClusterResults[[i]]) <- paste(rownames(deClusterResults[[i]]), cluster_ids[i], sep="-")
    }

    return(Reduce(rbind, deClusterResults))
}
