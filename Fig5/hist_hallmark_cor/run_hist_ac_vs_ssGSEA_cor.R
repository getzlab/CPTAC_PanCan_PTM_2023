library(biomaRt)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)



#load dendro mappings
mappings <-read.delim("../../data/full_dendro_sig_mappings.tsv")

colnames(mappings)<-c("Sample","d1","d3","d17","d22","d39","d54","d56","d5","d11","d29","d40","d24","d36","d13","d38","d19","d2","d6","d9","d10","d14","d18","d28","d20","d31","d25","d43","d23","d37","d33","d4","d7","d15","d21","d32","d8","d16","d26","d34","d12","max_id","max_norm","cohort")

dendro <-c("d1","d3","d17","d22","d39","d54","d56","d5","d11","d29","d40","d24","d36","d13","d38","d19","d2","d6","d9","d10","d14","d18","d28","d20","d31","d25","d43","d23","d37","d33","d4","d7","d15","d21","d32","d8","d16","d26","d34","d12")

for(g in dendro) {
    mappings[,g][mappings[,g]!=""] <- 1
}

d1_samples <- mappings %>% filter(d1==1) %>% pull(Sample)

##select histone sites
histone_genes <- read.delim('../../data/Histone_genes_v1.tsv')

acetyl <- read.delim("../../data/acetylome_imputed.tsv")


mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

refseq_ids <- unique(acetyl$short_refseq_id)

IDs <- getBM(attributes = c("hgnc_id","refseq_peptide"),
             filters = "refseq_peptide", values = refseq_ids,
             mart = mart)

colnames(IDs)[2] <- 'short_refseq_id'

acetyl<-dplyr::left_join(acetyl,IDs,by='short_refseq_id')


hist_acetyl_full <- acetyl %>%
    merge(., histone_genes, by='hgnc_id') %>%
    filter(grepl('histones', gene_type))



##load ssGSEA data

ssgsea_df=read.delim('../../data/rna_z-combined.gct', skip=2)
common_samples=intersect(colnames(hist_acetyl_full), colnames(ssgsea_df))

hist_acetyl <- hist_acetyl_full[,c('X', common_samples)] %>%
    column_to_rownames('X') %>% t %>%
    as.data.frame %>% 
    rownames_to_column('Sample')


ssgsea_mat <- ssgsea_df[,c('id', common_samples)] %>% column_to_rownames('id') %>% t %>% as.data.frame %>% 
    rownames_to_column('Sample') %>% .[,-c(2:4)]

hist_ssgsea_df <- merge(hist_acetyl, ssgsea_mat, by='Sample')

hist_sites <- colnames(hist_acetyl[,-1])
hallmark_paths <- colnames(ssgsea_mat[,-1])

for (i in 2:112) { hist_ssgsea_df[,i] <- as.numeric(hist_ssgsea_df[,i]) }

hist_path_cor_df <- data.frame()

for (site in hist_sites) {
    for (path in hallmark_paths) {
        cor_test <- cor.test(hist_ssgsea_df[,site], hist_ssgsea_df[,path], method='spearman')
        
        cor_test_df <- data.frame(hist_site=site, hallmark_path=path, 
                                  spearman_coef=as.numeric(cor_test$estimate), spearman_pval=as.numeric(cor_test$p.value))
        hist_path_cor_df <- rbind(hist_path_cor_df, cor_test_df)
        
    }
}

## for each histone acetylsite, adjust p-values to FDRs across 50 hallmark pathways

for (site in hist_sites) {
    hist_path_cor_df$fdr[hist_path_cor_df$hist_site==site] <- p.adjust(hist_path_cor_df$spearman_pval[hist_path_cor_df$hist_site==site], 
                                                                       method='fdr')
    
}

hist_path_annot_df <- hist_path_cor_df
colnames(hist_path_annot_df)[1] <- 'X'
hist_path_annot_df <- merge(hist_path_annot_df, hist_acetyl_full[,c('X', 'hgnc_id', 'refseq_id')], by='X')

#map to correct gene symbols using varmap
varmap_v4 <- read.delim('../../data/var_map_full_v4.tsv')

colnames(varmap_v4)[5] <- 'refseq_id'
refseq_to_gene_df <- distinct(varmap_v4[,c('refseq_id', 'geneSymbol')])

hist_path_annot_df <- merge(hist_path_annot_df, refseq_to_gene_df, by='refseq_id', all.x=T)

write.table(x=hist_path_annot_df, 
            file='hist_impute_vs_hallmark_pathway_cor_df_060622.tsv',
            quote=F, row.names=F, sep='\t')

