library(dplyr)
library(circlize)
library(RColorBrewer)

setwd('~/Box/Ding_Lab/Projects_Current/CPTAC/CPTACIII/CPTAC3_pancancer_PTM/Analysis/202211_01_lasso_redo/')

x = readRDS('objects/dendro_tumor_raw_v4_impute_iter100_enzyme_acetyl_full_111322.Rds')

clusters <- names(x)

all_dendro_enzyme_acetyl <- data.frame()

for (clust in clusters) {
print(clust)
    
dendro_enzyme_acetyl <-x[[clust]] %>%
    as.data.frame %>%
    mutate(lasso_min_l_MSE=as.numeric(lasso_min_l_MSE)) %>%
    filter(Enzyme_gene!='(Intercept)') %>%
    group_by(Enzyme_gene, Substrate_site) %>%
    mutate(min_MSE=min(lasso_min_l_MSE)) %>%
    mutate(lasso_coef=as.numeric(lasso_coef)) %>%
    filter(Enzyme_gene!=Substrate_site,lasso_min_l_MSE==min_MSE) %>%
    mutate(substrate_gene=strsplit(Substrate_site, split='\\.NP') %>% lapply(., '[[', 1) %>% unlist) %>%
    mutate(substrate_acetylsite=ifelse(grepl('\\.K', Substrate_site),strsplit(Substrate_site, split='\\.K') %>% lapply(., '[[', 2) %>% unlist,'')) %>%
    mutate(substrate_acetylsite=ifelse(substrate_acetylsite!='', paste0('K', substrate_acetylsite), substrate_acetylsite)) %>%
    mutate(abbrev_substrate_id=ifelse(substrate_acetylsite!='', paste0(substrate_gene, '-', substrate_acetylsite), substrate_gene)) %>%
    group_by(abbrev_substrate_id, Enzyme_gene) %>%
    mutate(avg_lasso_coef=mean(lasso_coef, na.rm=T)) %>%
    as.data.frame %>%
    .[,c('Enzyme_gene', 'substrate_gene', 'substrate_acetylsite', 'abbrev_substrate_id', 'avg_lasso_coef')] %>%
    filter(Enzyme_gene!=substrate_gene) %>%
    distinct %>%
    mutate(cluster=clust)


all_dendro_enzyme_acetyl <- rbind(all_dendro_enzyme_acetyl, dendro_enzyme_acetyl)
}







## CREBBP, H4.16-K13kK17k (all_samples, C14)


test.master <-readRDS("objects/test.maste_raw_impute.Rds")
mappings <-read.delim("../../Resources/cluster_map_112022.tsv")



ac_plot_sites <- c('H2BC18.NP_001154806.1.K17kK21k')
pro_plot_ids <- c('CREBBP')

out_dir <- 'figures/'

for (group in c('all_samples','C22')) {
    
    if (group=='all_samples') {
        plot_ac_pro <- test.master[,c('cluster_id',pro_plot_ids, ac_plot_sites)] 
    } else {
        plot_ac_pro <- test.master[,c('cluster_id',pro_plot_ids, ac_plot_sites)] %>% filter(cluster_id==group)
    }
    
    p <- ggplot(plot_ac_pro, aes(x=CREBBP, y=H2BC18.NP_001154806.1.K17kK21k))+
        geom_point(size=1, alpha=0.5)+
        #scale_color_manual(values=c('True'='#311B92', 'False'='#F9A825'))+
        geom_smooth(method=lm, na.rm = TRUE, fullrange= T, se=F,
                    aes(group=1),colour="black", size=0.5)+
        labs(x='CBP', y='H2B-K12K13k')+
        theme_classic()
    
    pdf(paste0(out_dir, group, '_CBP_H2B_scatter.pdf'), width=4, height=3)
    print(p)
    dev.off()
    
}

### NCOA1, H2Bc3 K16K17k (C12)
ac_plot_sites <- c('H2BC3.NP_066406.1.K16kK17k')
pro_plot_ids <- c('NCOA1')

group='C12'
plot_ac_pro <- test.master[,c('cluster_id',pro_plot_ids, ac_plot_sites)] %>% filter(cluster_id==group)

p <- ggplot(plot_ac_pro, aes(x=NCOA1, y=H2BC3.NP_066406.1.K16kK17k))+
    geom_point(size=1, alpha=0.5)+
    #scale_color_manual(values=c('True'='#311B92', 'False'='#F9A825'))+
    geom_smooth(method=lm, na.rm = TRUE, fullrange= T, se=F,
                aes(group=1),colour="black", size=0.5)+
    labs(x='NCOA1', y='H2B-K16K17k')+
    theme_classic()

pdf(paste0(out_dir, group, '_NCOA1_H2B_scatter.pdf'), width=4, height=3)
print(p)
dev.off()


## BRD4, H2BC18-K17kK21k (C7)
ac_plot_sites <- c('H2BC18.NP_001154806.1.K17kK21k')
pro_plot_ids <- c('BRD4')

group='C7'
plot_ac_pro <- test.master[,c('cluster_id',pro_plot_ids, ac_plot_sites)] %>% filter(cluster_id==group)

p <- ggplot(plot_ac_pro, aes(x=BRD4, y=H2BC18.NP_001154806.1.K17kK21k))+
    geom_point(size=1, alpha=0.5)+
    #scale_color_manual(values=c('True'='#311B92', 'False'='#F9A825'))+
    geom_smooth(method=lm, na.rm = TRUE, fullrange= T, se=F,
                aes(group=1),colour="black", size=0.5)+
    labs(x='BRD4', y='H2B-K17K21k')+
    theme_classic()

pdf(paste0(out_dir, group, '_BRD4_H2B_scatter.pdf'), width=4, height=3)
print(p)
dev.off()

## TAF1 H2B-K21 (C9)

ac_plot_sites <- c('H2BC3.NP_066406.1.K21k')
pro_plot_ids <- c('TAF1')

group='C9'
plot_ac_pro <- test.master[,c('cluster_id',pro_plot_ids, ac_plot_sites)] %>% filter(cluster_id==group)

p <- ggplot(plot_ac_pro, aes(x=TAF1, y=H2BC3.NP_066406.1.K21k))+
    geom_point(size=1, alpha=0.5)+
    #scale_color_manual(values=c('True'='#311B92', 'False'='#F9A825'))+
    geom_smooth(method=lm, na.rm = TRUE, fullrange= T, se=F,
                aes(group=1),colour="black", size=0.5)+
    labs(x='TAF1', y='H2B-K21k')+
    theme_classic()

pdf(paste0(out_dir, group, '_TAF1_H2B_scatter.pdf'), width=4, height=3)
print(p)
dev.off()




