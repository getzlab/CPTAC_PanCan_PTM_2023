## lasso regression for variable selection
## This script take acetyl/protein as input, 
## output the non-zero lasso regression coefficients between histone sites acetylation level and protein level
library (ggplot2)
library (permute)
library(data.table)
library(stringr)
library(glmnet)
library(corrplot)
library(gtools)
library(RColorBrewer)
#library(caret)
library(tidyverse)
library(impute)
library(dplyr)
library(biomaRt)
#library(SummarizedExperiment)



time=format(Sys.time(), "%Y%m%d") 
iter<-100 ## changing number of iterations to 100 
pval_enzyme_acetyl_raw<-NULL #List used to store all the boostrapping results from lasso regression
pval_enzyme_acetyl<-NULL #List used to store all non-zero lasso coefficients
rescoeff<-NULL #List used to store all linear regression coefficients

#add histone gene info, using Broad's varmap file now

varmap_v4 <- read.delim('../../../data/var_map_full_v4.tsv')
ac_site_to_gene_df <- varmap_v4 %>% filter(feature=='acetylome')

exclude_samples = c('Burdenko..1360_WNT', 'C3L.04350', 'C3L.05257', 'MB018_GR3', 'MB037_SHH', 'MB282_GR4', 'X03BR011', 'X1M6_WNT', 'X5M15_WNT')

#### check with Yizhe re: most recent version 
histone_genes <- read.delim('../../../data/Histone_genes_v1.tsv')
histone_genes <- histone_genes[,c('hgnc_id', 'gene_type', 'symbol')]
colnames(histone_genes)[3] <- 'geneSymbol'

ac_site_to_gene_df <- ac_site_to_gene_df %>%
  mutate(refseq_id=id %>% strsplit(., split='_') %>% lapply(., '[[', 2) %>% unlist %>% paste0('NP_',.)) %>%
  mutate(short_refseq_id =str_replace(refseq_id, '\\.([0-9]+)', ''))

ac_refseq_ids <- unique(ac_site_to_gene_df$short_refseq_id)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

ac_ids <- getBM(attributes = c("hgnc_id","refseq_peptide"),
                filters = "refseq_peptide", values = ac_refseq_ids,
                mart = mart)

colnames(ac_ids)[2] <- 'short_refseq_id'

ac_site_to_hgnc_df <- merge(ac_site_to_gene_df, ac_ids, by='short_refseq_id')

hist_ac_site_to_hgnc_df <- ac_site_to_hgnc_df %>%
  merge(., histone_genes[,-3], by='hgnc_id') %>%
  mutate(site=strsplit(id, split='_') %>% lapply(., '[[', 3) %>% unlist) %>%
  mutate(abbrev_site=paste0(geneSymbol, '-', refseq_id, '-', site))

#map acetyl to histone sites

mat_acetyl <- read.delim('../../../data/acetylome_X.tsv.gz')

colnames(mat_acetyl)[1] <- 'id'

mat_acetyl <- mat_acetyl %>%
  merge(., hist_ac_site_to_hgnc_df, by='id') %>%
  .[,c(597, 2:580)] %>%
  column_to_rownames('abbrev_site')

#map protein to histone genes/ids

mat_pro <- read.delim('../../../data/proteome_X.tsv.gz')

pro_refseq_ids <- unique(mat_pro$X)

pro_refseq_id_df <- data.frame(refseq_id=pro_refseq_ids)

pro_refseq_id_df$short_id <- str_replace(pro_refseq_id_df$refseq_id, '\\.([0-9]+)', '')

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

pro_ids <- getBM(attributes = c("hgnc_id","refseq_peptide"),
                 filters = "refseq_peptide", values = unique(pro_refseq_id_df$short_id),
                 mart = mart)

colnames(pro_ids)[2] <- 'short_id'

hist_pro_site_to_hgnc_df <-pro_ids %>%
  merge(., histone_genes, by='hgnc_id') %>%
  merge(., pro_refseq_id_df, by='short_id')

colnames(mat_pro)[1] <- 'refseq_id'

histone_pro <- mat_pro %>% 
  merge(., hist_pro_site_to_hgnc_df, by='refseq_id') %>%
  .[,c(1124,2:1120)]

histone_pro <- aggregate(histone_pro[,-1], list(gene=histone_pro[,1]), FUN=mean, na.rm=T)




# Read in acetylation data at Pancan Level
mat1 <-mat_acetyl

# Only keep sites with at least 80% of non-missing data
# You can change this cutoff
mat1 <-mat1[which(rowMeans(!is.na(mat1)) > 0.8), ]
# Impute NA with knn methods
library(impute)
mat2 <-impute.knn(as.matrix(mat1) ,k = 10, rowmax = 0.5, colmax = 0.8, maxp = 1500, rng.seed=362436069)


mat3 <-mat2
mat1<-mat3
rm(mat2)
rm(mat3)
saveRDS(mat1,"objects/mat1_raw_imputed.Rds")

mat1 <-readRDS("objects/mat1_raw_imputed.Rds")
mat1 <- mat1$data
# Read in protein data at Pancan Level
mat0 <-histone_pro %>% column_to_rownames('gene')
mat0 <- mat0[,colnames(mat1)]

# Combine protein/acetylation datasets
test.mat <- cbind(t(mat0),
                  t(mat1))
#colnames(test.mat)
#colnames(as_tibble(test.mat))

test.mat <- data.frame(case_id = row.names(test.mat), test.mat)
test.mat <- test.mat %>% filter(!(case_id %in% exclude_samples))

mappings <-read.delim("../../../data/cluster_map.tsv")


# Test if the case_id in test.mat and mappings are the same
colnames(mappings)[1] <- 'case_id'

test.master <- dplyr::left_join(test.mat,mappings,by="case_id")
saveRDS(test.master,"objects/test.maste_raw_impute.Rds")

# Read in test.master dataset
test.master <-readRDS("objects/test.maste_raw_impute.Rds")
## Prepare the matrix for lasso regression
test.mat <-NULL

## Lasso regression & Linear regression
## Select the proteins that has a sig impact on histone acetylation, the coefficient of non-sig protein will be shrinked to zero 
clusters <- unique(mappings$cluster_id)
clusters <- c(clusters, 'all_samples')

for (j in clusters) {
  
  exclude_cols <- c('case_id', 'cluster_id', 'max_id', 'max_norm', 'cohort')
  
  exclude_col_idx <- which(colnames(test.master) %in% exclude_cols)
  
  if (j=='all_samples') {
    test.mat <- test.master %>% filter(case_id %in% mappings$case_id) %>% .[,-exclude_col_idx]
  } else {
    test.mat <- test.master %>% filter(cluster_id==j) %>% .[,-exclude_col_idx]
  }
  rownames(test.mat) <- NULL
  
  site_cols <- which(colSums(is.na(test.mat))==0)
  test.mat <- test.mat[,site_cols]
  
  if (dim(test.mat)[1]!=0) {
    big.df.pval <- NULL
    res.pval <- NULL
    for (i in 1:iter) {
      df.pval <- NULL
      print(i)
      set.seed(i)
      # 80% training set, 20% testing set
      train = sample(1:nrow(test.mat), 4*nrow(test.mat)/5)
      test = (-train)
      
      for (site in colnames(test.mat)) {
        new <-test.mat
        new<-new[,!colnames(new)%in% colnames(data.frame(t(mat1)))]
        x <- data.matrix(new)
        ix <- which(colnames(test.mat) %in% site)
        y <- test.mat[,ix]
        ytest = y[test]
        # cross-validation
        cv.out <- cv.glmnet(x[train,], y[train], alpha=1, nfolds = 10, type.measure = 'mse')
        #min value of lambda
        lambda_min <- cv.out$lambda.min
        #best value of lambda
        lambda_1se <- cv.out$lambda.1se
        #regression coefficients
        lasso.coef <- coef(cv.out, s=lambda_1se)
        lasso.pred.se <- predict(cv.out, newx = x[test,], s=lambda_1se)
        lasso.pred.min <- predict(cv.out, newx = x[test,], s=lambda_min)
        df.pval <- rbind(df.pval, cbind(Enzyme_gene = rownames(lasso.coef),
                                        Substrate_site = site,
                                        N = nrow(x),
                                        best_lambda = lambda_1se,
                                        min_lambda = lambda_min,
                                        lasso_coef = lasso.coef[,1],
                                        lasso_best_l_MSE = mean((lasso.pred.se-ytest)^2),
                                        lasso_min_l_MSE = mean((lasso.pred.min-ytest)^2)))
        
      }
      big.df.pval <- rbind (big.df.pval, cbind(df.pval, i = i))
    }
    pval_enzyme_acetyl_raw[[j]]<-big.df.pval
    
    big.df.pval <- data.frame(big.df.pval)
    big.df.pval$enzyme_site <- paste(big.df.pval$Enzyme_gene, big.df.pval$Substrate_site, sep = ':')
    big.df.pval[,3:8] <- apply(big.df.pval[,3:8],MARGIN = c(1,2), as.numeric)
    big.df.pval <-big.df.pval[which(big.df.pval$lasso_coef != 0), ]
    
    pval_enzyme_acetyl[[j]]<-big.df.pval
    
    # This section calcultes the linear regression coefficient between sites~proteins at each dendro group
    for (site in colnames(test.mat)) {
      for ( i in colnames(test.mat)){
        form=paste0(site,'~',i)
        model <-lm(formula(form),data=test.mat)
        beta=model$coefficients[2]
        
        res.pval <- rbind(res.pval, cbind(Enzyme_gene = i,
                                          Substrate_site = site,
                                          beta=beta))
      }
    }
    rescoeff[[j]]=res.pval
  }
}

saveRDS( pval_enzyme_acetyl,"objects/dendro_tumor_raw_v4_impute_iter100_pval_enzyme_acetyl_111322.Rds")
saveRDS( pval_enzyme_acetyl_raw,"objects/dendro_tumor_raw_v4_impute_iter100_enzyme_acetyl_full_111322.Rds")
saveRDS( rescoeff,"objects/dendro_tumor_raw_v4_impute_iter100_rescoeff_111322.Rds")

print('done')






