
raw_acetyl=read.delim('../../data/acetylome_raw_res_X.tsv')
raw_phospho=read.delim('../../data//phosphoproteome_raw_res_X.tsv')
site_col_name='id'

raw_acetyl <- raw_acetyl %>% column_to_rownames(site_col_name) %>%
  t %>% as.data.frame %>%
  rownames_to_column('Sample')

raw_phospho <- raw_phospho %>% column_to_rownames(site_col_name) %>%
  t %>% as.data.frame %>%
  rownames_to_column('Sample')

raw_ac_pho <- merge(raw_acetyl, raw_phospho, by='Sample')

#format acetyl sites
all_acetyl_site_df <- data.frame(ac_site_id=unique(colnames(raw_acetyl)[-1]))
all_acetyl_site_df <- all_acetyl_site_df %>%
  mutate(ac_sites=ac_site_id %>% strsplit(., split='_') %>% lapply(., '[[', 3) %>% unlist) %>%
  mutate(ac_sep_sites=strsplit(ac_sites, 'k') ) %>%
  tidyr::unnest(ac_sep_sites) %>%
  mutate(ac_loc=tidyr::extract_numeric(ac_sep_sites)) %>%
  mutate(refseq_id_prefix = strsplit(ac_site_id, split='_') %>% lapply(., '[[', 1) %>% unlist) %>%
  mutate(refseq_id_num=strsplit(ac_site_id, split='_') %>% lapply(., '[[', 2) %>% unlist) %>%
  mutate(refseq_id=paste0(refseq_id_prefix, '_', refseq_id_num))

#format phospho sites
all_phospho_site_df <- data.frame(pho_site_id=unique(colnames(raw_phospho)[-1]))
all_phospho_site_df <- all_phospho_site_df %>%
  mutate(pho_sites=pho_site_id %>% strsplit(., split='_') %>% lapply(., '[[', 3) %>% unlist) %>%
  mutate(pho_sep_sites=strsplit(pho_sites, 's|t|y')) %>%
  tidyr::unnest(pho_sep_sites) %>%
  mutate(pho_loc=tidyr::extract_numeric(pho_sep_sites)) %>%
  mutate(refseq_id_prefix = strsplit(pho_site_id, split='_') %>% lapply(., '[[', 1) %>% unlist) %>%
  mutate(refseq_id_num=strsplit(pho_site_id, split='_') %>% lapply(., '[[', 2) %>% unlist) %>%
  mutate(refseq_id=paste0(refseq_id_prefix, '_', refseq_id_num))


all_adj_ac_pho_site_df <- data.frame()

for (i in 1:nrow(all_acetyl_site_df)) {
  ac_site_df <- all_acetyl_site_df[i,]
  pho_site_df <- all_phospho_site_df %>% 
    filter(refseq_id==ac_site_df$refseq_id, abs(pho_loc-ac_site_df$ac_loc) <= 5)
  
  if (nrow(pho_site_df) > 0) {
    merge_site_df <- merge(ac_site_df, pho_site_df[,-c(5:6)], by='refseq_id', all.y=T)
  } else {
    na_pho_site_df <- data.frame(pho_site_id=NA, pho_sites=NA, pho_sep_sites=NA, pho_loc=NA)
    ac_site_df <- ac_site_df[,c(7,1:6)]
    merge_site_df <- cbind(ac_site_df, na_pho_site_df)
  }
  
  all_adj_ac_pho_site_df <- rbind(all_adj_ac_pho_site_df, merge_site_df)
  print(i/nrow(all_acetyl_site_df))
}

write.table(x=all_adj_ac_pho_site_df,
            file=paste0('tables/all_adj_ac_pho_sites_', group, '_10222021.tsv'),
            quote=F, sep='\t', row.names=F)
print('wrote site file')