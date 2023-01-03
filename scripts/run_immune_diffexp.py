from diffexp_wrapper import run_differential_expression, load_de, load_mut_de, gen_phospho_ptmgsea, load_immune_de
import glob
import pandas as pd
import os

SCRIPTS_DIR = "."

# ---------------------------------------------------------
# Run Diffexp for RNA, Protein, Phospho per Dendro Group w/ Missing
# ---------------------------------------------------------
mappings_file = "../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune/mappings.tsv"
mdf = pd.read_csv(mappings_file, sep='\t', index_col=0)

for col in mdf.columns:
    print("Running for group {}...".format(col))
    formats = dict()
    formats['input'] = mappings_file
    formats['clust'] = col
    formats['proteome']  = "../data/processed/061721_genecentric/raw/proteome_X.tsv.gz"
    formats['phosphoproteome']  = "../data/processed/061721/raw_res/phosphoproteome_raw_res_X.tsv"
    formats['transcriptome'] = "../data/processed/061721/rna/tumor_rna_counts_pc_X.tsv"
    formats['acetylome']  = "../data/processed/061721/raw_res/acetylome_raw_res_X.tsv"
    formats['prot_maps'] = "../data/processed/061721/var_map_full.tsv"
    formats['covar'] = "../analysis/signatures/061721_imputed_res_reg/mappings.tsv"
    formats['out_limma'] = "../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune/immune_{}_cohort_cov".format(col)
    formats["minObs"] = 10

    if col =="BestCall":
        formats["subset"] = "IFN-Gamma Inflammatory Lymphocyte_Depleted TGF-Beta Wound_Healing"
    run_differential_expression(formats, scripts_dir=SCRIPTS_DIR)

full_dendro_de_cohort = pd.concat([load_immune_de(x) for x in glob.glob("../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune/*cohort_cov/*full*")])

# Added for proteome gene-centric
full_dendro_de_cohort.loc[full_dendro_de_cohort['feature'].isna(),'gene_name'] = full_dendro_de_cohort.loc[full_dendro_de_cohort['feature'].isna(),'index']
full_dendro_de_cohort.loc[full_dendro_de_cohort['feature'].isna(),'feature'] = 'proteome'

# Save
full_dendro_de_cohort.to_csv("../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune/full_de_cohort_cov.tsv",sep='\t')

pmap_df = pd.read_csv("../data/processed/061721/var_map_full.tsv", sep='\t', index_col=0)
de_df = pd.read_csv("../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune/full_de_cohort_cov.tsv", sep='\t', index_col=0).set_index('index')
gen_phospho_ptmgsea(
    de_df,
    pmap_df,
    "../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune/de_phosph_ptmgsea_cohort_cov.gct"
)

# ---------------------------------------------------------
# Not Residualized
# ---------------------------------------------------------
os.makedirs("../analysis/Fig_immuno_metabolism/diffexp/061721_raw_immune/", exist_ok=True)

mappings_file = "../analysis/Fig_immuno_metabolism/diffexp/061721_raw_res_immune/mappings.tsv"
mdf = pd.read_csv(mappings_file, sep='\t', index_col=0)

for col in mdf.columns:
    print("Running for group {}...".format(col))
    formats = dict()
    formats['input'] = mappings_file
    formats['clust'] = col
    formats['phosphoproteome']  = "../data/processed/061721/raw/phosphoproteome_X.tsv.gz"
    formats['acetylome']  = "../data/processed/061721/raw/acetylome_X.tsv.gz"
    formats['prot_maps'] = "../data/processed/061721/var_map_full.tsv"
    formats['covar'] = "../analysis/signatures/061721_imputed_res_reg/mappings.tsv"
    formats['out_limma'] = "../analysis/Fig_immuno_metabolism/diffexp/061721_raw_immune/immune_{}_cohort_cov".format(col)
    formats["minObs"] = 10

    if col =="BestCall":
        formats["subset"] = "IFN-Gamma Inflammatory Lymphocyte_Depleted TGF-Beta Wound_Healing"

    run_differential_expression(formats, scripts_dir=SCRIPTS_DIR)

full_dendro_de_cohort = pd.concat([load_immune_de(x) for x in glob.glob("../analysis/Fig_immuno_metabolism/diffexp/061721_raw_immune/*cohort_cov/*full*")])

# Save
full_dendro_de_cohort.to_csv("../analysis/Fig_immuno_metabolism/diffexp/061721_raw_immune/full_de_cohort_cov.tsv",sep='\t')

pmap_df = pd.read_csv("../data/processed/061721/var_map_full.tsv", sep='\t', index_col=0)
de_df = pd.read_csv("../analysis/Fig_immuno_metabolism/diffexp/061721_raw_immune/full_de_cohort_cov.tsv", sep='\t', index_col=0).set_index('index')
gen_phospho_ptmgsea(
    de_df,
    pmap_df,
    "../analysis/Fig_immuno_metabolism/diffexp/061721_raw_immune/de_phosph_ptmgsea_cohort_cov.gct"
)
