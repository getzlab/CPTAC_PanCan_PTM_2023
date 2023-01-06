
library(tidyverse)
library(dplyr)
library(biomaRt)


######## with protein-corrected ptm abundance

raw_acetyl=read.delim('../../data/acetylome_raw_res_X.tsv')
raw_phospho=read.delim('../../data/phosphoproteome_raw_res_X.tsv')

mappings <-read.delim("../../data/cluster_map.tsv")

colnames(mappings)[1] <- 'Sample'

clusters <- unique(mappings$cluster_id)

raw_acetyl <- raw_acetyl %>% column_to_rownames('id') %>%
    t %>% as.data.frame %>%
    rownames_to_column('Sample')

raw_phospho <- raw_phospho %>% column_to_rownames('id') %>%
    t %>% as.data.frame %>%
    rownames_to_column('Sample')

map_acetyl <- merge(raw_acetyl, mappings, by='Sample')

map_ac_pho <- merge(map_acetyl, raw_phospho, by='Sample')

exclude_samples = c('Burdenko..1360_WNT', 'C3L.04350', 'C3L.05257', 'MB018_GR3', 'MB037_SHH', 'MB282_GR4', 'X03BR011', 'X1M6_WNT', 'X5M15_WNT')

map_ac_pho = map_ac_pho %>% filter(!(Sample %in% exclude_samples))
print(dim(map_ac_pho))

ac_pho_site_df <- read.delim('../../data/all_adj_ac_pho_sites_ptm_protein_corrected_10222021.tsv')

clusters <- c(clusters, 'all_samples')

all_crosstalk_df <- data.frame()



for (i in c('all_samples',1:23)) {
    
    if (i!='all_samples') { 
        i=as.numeric(i) 
        clusts = paste0('C',i:(i+1))
    } else {
        clusts = i
    }
    
    
    clust_crosstalk_df <- data.frame()
    
    if (i=='all_samples') {
        samples <- mappings$Sample
    } else {
        samples <- mappings$Sample[mappings$cluster_id %in% clusts]
    }
    
    clust_ac_pho <- map_ac_pho %>% filter(Sample %in% samples)
    for (i in 1:nrow(ac_pho_site_df)) {
        df <- ac_pho_site_df[i,]
        
        df$cor <- NA
        df$p.value <- NA
        
        ac_site <- df$ac_site_id
        pho_site <- df$pho_site_id
        
        if (all(c(ac_site, pho_site) %in% colnames(clust_ac_pho))) {
            ac_pho_val_df <- clust_ac_pho[,c(ac_site, pho_site)]
            
            non_na_idx <- which(!is.na(ac_pho_val_df[,1]) & !is.na(ac_pho_val_df[,2]))
            if (length(non_na_idx) >= 3) {
                tst <- cor.test(ac_pho_val_df[,1], ac_pho_val_df[,2], use='p', method='spearman')
                df$cor <- tst$estimate; df$p.value <- tst$p.value
            }
        }
        
        print(paste(clusts, i/nrow(ac_pho_site_df), collapse=' '))
        
        clust_crosstalk_df <- rbind(clust_crosstalk_df, df)
    }
    
    
    clust_crosstalk_df$cluster <- paste(clusts, collapse='/')
    first_clust=clusts[1]
    write.table(x=clust_crosstalk_df, file=paste0('tables/', first_clust, '_adj_clust_crosstalk_cor_protein_corrected_122722.tsv'), quote=F, row.names=T, sep='\t')
    all_crosstalk_df <- rbind(all_crosstalk_df, clust_crosstalk_df)
    
    
}

write.table(x=all_crosstalk_df,
            file='tables/all_crosstalk_cor_adjacent_clusters_protein_corrected_12272022.tsv',
            quote=F, row.names=F, sep='\t')

