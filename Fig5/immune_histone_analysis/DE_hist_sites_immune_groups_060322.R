


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
    filter(grepl('histones', gene_type)) %>%
    .[,c(2:581)] %>%
    column_to_rownames('X') %>% t %>%
    as.data.frame %>%
    rownames_to_column('Sample')
    

#immune mappings

immune_map_df <- read.delim('../../data/mappings_immune.tsv')
colnames(immune_map_df)[1] <- 'Sample'

#Cool vs all others comparison

for (group in unique(immune_map_df$ImmuneClustC4)) {

    print(group)
    
    im_group_samples <- immune_map_df$Sample[immune_map_df$ImmuneClustC4==group]
    
    
    hist_im_acetyl <- hist_acetyl_full %>% filter(Sample %in% im_group_samples) %>% column_to_rownames('Sample')
    
    hist_im_not_acetyl <- hist_acetyl_full %>% filter(!(Sample %in% im_group_samples)) %>% column_to_rownames('Sample')
    
    #run Wilcoxon test
    hist_sites <- colnames(hist_acetyl_full)[-1]
    
    hist_group_vs_others_df <- data.frame(site=hist_sites, statistic=NA, p_val=NA, median=NA, median_not=NA, median_diff=NA)
    group_name = tolower(str_replace(group, 'Immune', ''))
    colnames(hist_group_vs_others_df)[4:5] <- paste0(colnames(hist_group_vs_others_df)[4:5], '_', group_name)
    
    for (i in 1:nrow(hist_group_vs_others_df)) {
        
        print(i)
        hist_site=hist_group_vs_others_df$site[i]
        site_ac_group <- hist_im_acetyl[,hist_site]
        site_ac_not_group <- hist_im_not_acetyl[,hist_site]
        
        tst <- wilcox.test(x=site_ac_group, y=site_ac_not_group, paired=F)
        
        hist_group_vs_others_df$statistic[i] <- tst$statistic
        hist_group_vs_others_df$p_val[i] <- tst$p.value
        hist_group_vs_others_df$median_group[i] <- median(site_ac_group, na.rm=T)
        hist_group_vs_others_df$median_not_group[i] <- median(site_ac_not_group, na.rm=T)
        
        hist_group_vs_others_df$median_diff[i] <- hist_group_vs_others_df$median_group[i] - hist_group_vs_others_df$median_not_group[i]
            
    }
    
    hist_group_vs_others_df$fdr <- p.adjust(hist_group_vs_others_df$p_val, method='fdr')
    
    #map to varmap
    varmap_v4 <- read.delim('../../data/var_map_full_v4.tsv')
    colnames(varmap_v4)[1] <- 'site'
    
    hist_group_vs_others_df <- merge(hist_group_vs_others_df, varmap_v4[,c('site', 'geneSymbol')], by='site')
    
    write.table(x=hist_group_vs_others_df,
                file='hist_DE_immune_', group_name,'_vs_others_062522.tsv',
                quote=F, row.names=F, sep='\t')
    
}


