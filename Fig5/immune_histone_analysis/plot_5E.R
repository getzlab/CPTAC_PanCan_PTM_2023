library(biomaRt)
library(dplyr)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)


out_dir = '~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Figures/figure_panels/panel_E/'
setwd(out_dir)

##select histone sites
histone_genes <- read.delim('~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Analysis/forAkshay/testset/histone_related_genes.tsv')
hgnc_data <- read.delim('~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Data_Freeze/Broad/v1/Analyses_files_061721/mappings/HGNC_database.txt')


acetylome_imputed <- read.delim("~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Data_Freeze/Broad/imputed_from_drive/acetylome_X.tsv.gz")


acetyl<-acetylome_imputed %>%
    mutate(refseq_prefix=strsplit(X, split='_') %>% lapply(., '[[', 1) %>% unlist) %>%
    mutate(refseq_num=strsplit(X, split='_') %>% lapply(., '[[', 2) %>% unlist) %>%
    mutate(short_refseq_num=strsplit(refseq_num, split='\\.') %>% lapply(., '[[', 1) %>% unlist) %>%
    mutate(short_refseq_id=paste(refseq_prefix, short_refseq_num, sep='_')) %>%
    mutate(refseq_id=paste(refseq_prefix, refseq_num, sep='_'))

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

ssgsea_df=read.delim('~/Downloads/rna_z-combined.gct', skip=2)
common_samples=intersect(colnames(hist_acetyl_full), colnames(ssgsea_df))

hist_acetyl <- hist_acetyl_full[,c('X', common_samples)] %>%
    column_to_rownames('X') %>% t %>%
    as.data.frame %>% 
    rownames_to_column('Sample')


ssgsea_mat <- ssgsea_df[,c('id', common_samples)] %>% column_to_rownames('id') %>% t %>% as.data.frame %>% 
    rownames_to_column('Sample') %>% .[,-c(2:4)]

ssgsea_hist_mat <- merge(hist_acetyl, ssgsea_mat, by='Sample')

hist_sites <- colnames(hist_acetyl[,-1])
hallmark_paths <- colnames(ssgsea_mat[,-1])


for (i in 2:112) { ssgsea_hist_mat[,i] <- as.numeric(ssgsea_hist_mat[,i]) }

#map dendro groups
mappings <-read.delim("~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Analysis/forAkshay/ResultsFromGadGroup/full_dendro_sig_mappings.tsv")

colnames(mappings)<-c("Sample","d1","d3","d17","d22","d39","d54","d56","d5","d11","d29","d40","d24","d36","d13","d38","d19","d2","d6","d9","d10","d14","d18","d28","d20","d31","d25","d43","d23","d37","d33","d4","d7","d15","d21","d32","d8","d16","d26","d34","d12","max_id","max_norm","cohort")

dendro <-c("d1","d3","d17","d22","d39","d54","d56","d5","d11","d29","d40","d24","d36","d13","d38","d19","d2","d6","d9","d10","d14","d18","d28","d20","d31","d25","d43","d23","d37","d33","d4","d7","d15","d21","d32","d8","d16","d26","d34","d12")

for(g in dendro) {
    mappings[,g][mappings[,g]!=""] <- 1
}

ssgsea_hist_mat <- merge(ssgsea_hist_mat, mappings, by='Sample')

ssgsea_hist_mat$dendro[ssgsea_hist_mat$d4==1]<-"dendro 4"
ssgsea_hist_mat$dendro[ssgsea_hist_mat$d5==1]<-"dendro 5"
ssgsea_hist_mat$dendro[ssgsea_hist_mat$d6==1]<-"dendro 6"
ssgsea_hist_mat$dendro[ssgsea_hist_mat$d17==1]<-"dendro 17"
ssgsea_hist_mat$dendro[is.na(ssgsea_hist_mat$dendro)]<-"NA"

# map gene symbols
varmap_v4 <- read.delim('~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Data_Freeze/Broad/v1/Analyses_files_061721/processed_Data/061721/var_map_full_v4.tsv')


select_pths <- c('HALLMARK_GLYCOLYSIS', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_E2F_TARGETS', 
                 'HALLMARK_MYC_TARGETS_V2', 'HALLMARK_G2M_CHECKPOINT')


plot_hist_ssgsea <- ssgsea_hist_mat[,c('Sample', 'cohort', 'NP_003520.1_K28kK37k_2_2_28_37', 'NP_003520.1_K37k_1_1_37_37', select_pths)]
colnames(plot_hist_ssgsea)[3:4] <- c('H3C12-K28K37', 'H3C12-K37')

##add immune mappings
immune_map_df <- read.delim('~/Downloads/mappings_labeled (1).tsv')
colnames(immune_map_df)[1] <- 'Sample'

plot_hist_ssgsea=merge(plot_hist_ssgsea, immune_map_df[,c('Sample', 'ImmuneClustC4')], by='Sample')

col_cutoff=plot_hist_ssgsea[,-c(1:2, ncol(plot_hist_ssgsea))] %>% scale %>% t %>% abs %>% max(., na.rm=T)

col_fun = colorRamp2(c(-3, 0, 3), rev(brewer.pal(3,"RdBu")))

plot_hist_ssgsea$cohort[plot_hist_ssgsea$cohort=='MEDUL'] <- 'MB'
plot_hist_ssgsea$ImmuneClustC4 <- factor(plot_hist_ssgsea$ImmuneClustC4, levels=c('ImmuneCold', 'ImmuneCool', 'ImmuneWarm', 'ImmuneHot'))

h1=Heatmap(plot_hist_ssgsea[,c('H3C12-K28K37', 'H3C12-K37', select_pths)] %>% scale %>% t,
           col=col_fun,
           name='Z score',
           row_names_gp=gpar(fontsize=9),
           column_split=plot_hist_ssgsea$ImmuneClustC4,
           cluster_column_slices = F,
           bottom_annotation=columnAnnotation(Cohort=plot_hist_ssgsea$cohort,
                                              Immune=plot_hist_ssgsea$ImmuneClustC4,
                                              col=list(Cohort=c("GBM"="#52b788", "MB"="#193e2e", "LSCC"="#91bdff", "LUAD"="#1a759f", "UCEC"="#5a189a", "BRCA"="#cd6090"),
                                                       Immune=c('ImmuneCool'='#73ABC9', 'ImmuneWarm'='#F6A43E', 'ImmuneCold'='#8A74C0', 'ImmuneHot'='#B94A47'))))
pdf(paste0(out_dir, 'hist_acetyl_imputed_selected_ssGSEA_pathways_by_immune_060222.pdf'), width=18, height=3)
h1
dev.off()

################## sites decreased in immune cool (raw acetyl) #########################

##check histone sites decreased in immune

immune_de_df <- read.delim('~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Data_Freeze/Broad/immune_DE_results/061721_raw_res_immune_ciber_subtypes_named.tsv')
rownames(immune_de_df) <- NULL

#map refseq accessions to HGNC IDs
immune_de_df$short_refseq_id <- immune_de_df$accession_number %>% str_replace('\\.[0-9]+', '')

refseq_to_hgnc_df <- getBM(attributes = c("hgnc_id","refseq_peptide"),
                           filters = "refseq_peptide", values = unique(immune_de_df$short_refseq_id),
                           mart = mart)

colnames(refseq_to_hgnc_df)[2] <- 'short_refseq_id'

immune_annot_de_df <- merge(immune_de_df, refseq_to_hgnc_df, by='short_refseq_id', all.x=T)
hist_immune_annot_de_df <- merge(immune_annot_de_df, histone_genes, by='hgnc_id') 

im_cool_sites <- hist_immune_annot_de_df %>%
    filter(grepl('histones', gene_type), feature=='acetylome', qval < 0.1, id=='Cool', logFC < 0) %>%
    pull(index)

raw_acetyl <- read.delim('~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Data_Freeze/Broad/v1/Analyses_files_061721/processed_Data/061721/raw/acetylome_X.tsv.gz')

immune_map_df <- read.delim('~/Downloads/mappings_labeled (1).tsv')
colnames(immune_map_df)[1] <- 'Sample'

im_cool_site_acetyl <- raw_acetyl %>%
    filter(X %in% im_cool_sites) %>%
    mutate(id=X) %>%
    merge(., varmap_v4[,c('id', 'geneSymbol')], all.x=T) %>%
    mutate(site=strsplit(id, split='_') %>% lapply(., '[[', 3) %>% unlist) %>%
    mutate(short_site=paste(geneSymbol, site, sep='-')) %>%
    .[,c(584, 3:581)] %>%
    column_to_rownames('short_site') %>%
    t %>% as.data.frame %>%
    rownames_to_column('Sample') %>%
    merge(., immune_map_df[,c('Sample', 'ImmuneClustC4')], by='Sample') %>%
    merge(., mappings[,c('Sample', 'cohort')], by='Sample') %>%
    mutate(cohort=ifelse(cohort=='MEDUL', 'MB', cohort))

h2=Heatmap(im_cool_site_acetyl %>% .[,-c(1, ncol(.)-1, ncol(.))] %>% scale %>% t,
              col=col_fun,
              name='Z score',
              row_names_gp=gpar(fontsize=9),
              column_split=im_cool_site_acetyl$ImmuneClustC4,
              cluster_column_slices = F, 
              bottom_annotation=columnAnnotation(Cohort=im_cool_site_acetyl$cohort,
                                                 Immune=im_cool_site_acetyl$ImmuneClustC4,
                                                 col=list(Cohort=c("GBM"="#52b788", "MB"="#193e2e", "LSCC"="#91bdff", "LUAD"="#1a759f", "UCEC"="#5a189a", "BRCA"="#cd6090"),
                                                          Immune=c('ImmuneCool'='#73ABC9', 'ImmuneWarm'='#F6A43E', 'ImmuneCold'='#8A74C0', 'ImmuneHot'='#B94A47'))))

pdf(paste0(out_dir, 'hist_acetyl_raw_immune_cool_060222.pdf'), width=18, height=3)
h2
dev.off()

### combine two heatmaps

# same as h1 without annotations
h3 = Heatmap(plot_hist_ssgsea[,c('H3C12-K28K37', 'H3C12-K37', select_pths)] %>% scale %>% t,
                col=col_fun,
                name='Z score',
                row_names_gp=gpar(fontsize=9),
                column_split=plot_hist_ssgsea$ImmuneClustC4,
                cluster_column_slices = F)
    
pdf(paste0(out_dir, 'hist_acetyl_impute_ssGSEA_raw_immune_cool_060222.pdf'), width=18, height=5)
h3 %v% h2
dev.off()


################## sites decreased in immune cool (imputed acetyl) #########################

im_cool_impute_de_df <- read.delim('~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Tables/supp_tables/hist_DE_immune_cool_vs_others_060322.tsv')
im_cool_sites <- im_cool_impute_de_df %>% filter(fdr < 0.1, median_diff < 0) %>% pull(site)

im_cool_short_sites <- im_cool_impute_de_df %>%
    column_to_rownames('site') %>%
    .[im_cool_sites,] %>%
    rownames_to_column('site') %>% 
    mutate(loc=strsplit(site, split='_') %>% lapply(., '[[', 3) %>% unlist) %>%
    mutate(short_site=paste(geneSymbol, loc, sep='-')) %>%
    pull(short_site)

plot_hist_ssgsea <- ssgsea_hist_mat[,c('Sample', 'cohort', im_cool_sites)]
colnames(plot_hist_ssgsea)[-c(1:2)] <- im_cool_short_sites

##add immune mappings
immune_map_df <- read.delim('~/Downloads/mappings_labeled (1).tsv')
colnames(immune_map_df)[1] <- 'Sample'

plot_hist_ssgsea=merge(plot_hist_ssgsea, immune_map_df[,c('Sample', 'ImmuneClustC4')], by='Sample')

col_cutoff=plot_hist_ssgsea[,-c(1:2, ncol(plot_hist_ssgsea))] %>% scale %>% t %>% abs %>% max(., na.rm=T)

col_fun = colorRamp2(c(-3, 0, 3), rev(brewer.pal(3,"RdBu")))

plot_hist_ssgsea$cohort[plot_hist_ssgsea$cohort=='MEDUL'] <- 'MB'
plot_hist_ssgsea$ImmuneClustC4 <- factor(plot_hist_ssgsea$ImmuneClustC4, levels=c('ImmuneCold', 'ImmuneCool', 'ImmuneWarm', 'ImmuneHot'))

h4=Heatmap(plot_hist_ssgsea[,im_cool_short_sites] %>% scale %>% t,
           col=col_fun,
           name='Z score',
           row_names_gp=gpar(fontsize=9),
           column_split=plot_hist_ssgsea$ImmuneClustC4,
           cluster_column_slices = F,
           bottom_annotation=columnAnnotation(Cohort=plot_hist_ssgsea$cohort,
                                              Immune=plot_hist_ssgsea$ImmuneClustC4,
                                              col=list(Cohort=c("GBM"="#52b788", "MB"="#193e2e", "LSCC"="#91bdff", "LUAD"="#1a759f", "UCEC"="#5a189a", "BRCA"="#cd6090"),
                                                       Immune=c('ImmuneCool'='#73ABC9', 'ImmuneWarm'='#F6A43E', 'ImmuneCold'='#8A74C0', 'ImmuneHot'='#B94A47'))))

pdf(paste0(out_dir, 'hist_acetyl_impute_immune_cool_060322.pdf'), width=18, height=5)
h4
dev.off()

pdf(paste0(out_dir, 'hist_acetyl_impute_ssGSEA_immune_cool_060322.pdf'), width=18, height=5)
h3 %v% h4
dev.off()
