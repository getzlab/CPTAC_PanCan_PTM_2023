import numpy as np
import os
import pandas as pd
import glob
import sys
import matplotlib.pyplot as plt
import argparse

sys.path.append("/home/yakiyama/CPTAC_PanCan_2021/funcs")

import plotting as pl
import proteomics as prot

def main():
    parser = argparse.ArgumentParser(
        description='Process differential expression results.'
    )
    parser.add_argument(
        '-i', '--input_dir',
        required=True,
        help='<Required> Directory with differential expression outputs.'
    )
    parser.add_argument(
        '-f', '--feature_maps',
        required=True,
        help='<Required> File mapping features to genes.'
    )
    parser.add_argument(\
        '-o', '--out_dir',
        help='<Required> Output directory for plots.',
        required=True,
        type=str
    )

    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    map_df = pd.read_csv(args.feature_maps, sep='\t', index_col=0)

    # -----------------------
    # Aggregate DiffExp
    # -----------------------
    de_df = list()

    for de_file in glob.glob(os.path.join(args.input_dir, "*_de.tsv")):
        print("   * loading {}".format(de_file))
        _df = pd.read_csv(de_file, sep='\t', index_col=0).set_index("index")
        _df = _df.join(map_df[['id.description','variableSites','accession_number']])
        _df['feature'] = os.path.basename(de_file).split('_de')[0]
        _df['gsea_rank'] = _df['logFC'] * -np.log10(_df['adj.P.Val'])
        _df['gsea_rank_p'] = _df['logFC'] * -np.log10(_df['P.Value'])
        _df = _df.join(map_df[['causalpath_adjusted_id']])
        de_df.append(_df)

    de_df = pd.concat(de_df)
    de_df.to_csv(os.path.join(args.input_dir, "full_diffexp_results.tsv"), sep='\t')

    # # -----------------------
    # # Volcano plots
    # # -----------------------
    # de_df_prot = de_df[de_df['feature']!='Transcriptome']
    # de_df_rna = de_df[de_df['feature']=='Transcriptome']
    #
    # def _order(x):
    #     if x=='Acetylome':
    #         return 1
    #     elif x=='Proteome':
    #         return 2
    #     elif x=='Phosphoproteome':
    #         return 3
    #
    # de_df_prot["_order"] = de_df_prot['feature'].apply(_order).values
    #
    # # Hi-res plots with Genes
    # print("   * plotting hi-res protein/ptm volcanos")
    # os.makedirs(os.path.join(args.out_dir, "protein_ptm"), exist_ok=True)
    #
    # for idx in pd.unique(de_df_prot['id']):
    #     pl.plot_volcano(
    #         de_df_prot[de_df_prot['id']==idx].sort_values('_order', ascending=False),
    #         yax='adj.P.Val',
    #         xax='logFC',
    #         gene_id='gene_name',
    #         ptm_id='feature',
    #         thresh=2,
    #         shuffle=False,
    #         ylabel='-log10 (adj. p-val)',
    #         xlabel='log2 FC',
    #         gene_fontsize=5,
    #         only_plot_gene_symbol=False,
    #         label_percentile=99.925,
    #         c1='lightgrey'
    #     )
    #     plt.title("Cluster {}".format(idx), fontsize=16)
    #     plt.savefig(os.path.join(args.out_dir, "protein_ptm", "id_site_{}_de.pdf".format(idx)), dpi=200, bbox_inches='tight')
    #     plt.close()
    #
    #     os.makedirs(os.path.join(args.out_dir, "rna"), exist_ok=True)
    #     pl.plot_volcano(
    #         de_df_rna[de_df_rna['id']==idx],
    #             yax='adj.P.Val',
    #             xax='logFC',
    #             gene_id='gene_name',
    #             ylabel='-log10 (adj. p-val)',
    #             xlabel='log2 FC',
    #             ptm_id=None,
    #             arrow=False,
    #             gene_fontsize=5,
    #             thresh=2,
    #             label_percentile=99.925,
    #             c1='lightgrey',
    #             c2= pl.PTM_SCHEME['RNA']
    #         )
    #
    #     plt.title("Cluster {}".format(idx), fontsize=16)
    #     plt.savefig(os.path.join(args.out_dir, "rna", "id_rna_{}_de.pdf".format(idx)), dpi=200, bbox_inches='tight')
    #     plt.close()
    #
    # # -----------------------
    # # QQ Plots
    # # -----------------------
    # # # Volcano + QQ plot
    # print("   * plotting volcano + qq plots")
    # os.makedirs(os.path.join(args.out_dir, "vc_qq"), exist_ok=True)
    #
    # for feat in pd.unique(de_df['feature']):
    #     os.makedirs(os.path.join(args.out_dir, "vc_qq", feat), exist_ok=True)
    #
    #     pl.plot_de_qq_volcano(
    #         de_df[de_df['feature']==feat],
    #         groupby='id',
    #         xax='logFC',
    #         save_dir=os.path.join(args.out_dir, "vc_qq", feat),
    #         pval_idx='P.Value',
    #         yax='adj.P.Val',
    #         gene_id='gene_name',
    #         thresh=2,
    #         label_percentile=100,
    #         c1='lightgrey',
    #         c2=pl.PTM_SCHEME2[feat]
    #     )
    #     plt.close()

    # -----------------------
    # PTM Proportions
    # -----------------------
    print("   * plotting proportions")
    up_df, down_df = prot.get_ptm_de_proportions(
        de_df.drop(columns=['feature']),
        map_df['feature'],
        norm=False,
        lfc_thresh=0.5,
        qval_idx='adj.P.Val',
        lfc_idx='logFC',
        cluster_idx='id'
    )
    _ = pl.plot_ptm_de_proportions(up_df, down_df)
    plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, "ptm_proportions.pdf"), dpi=200, bbox_inches='tight')

if __name__ == "__main__":
	main()
