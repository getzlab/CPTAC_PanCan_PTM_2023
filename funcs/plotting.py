import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

import sys
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from typing import Union
import seaborn as sns
import signatureanalyzer as sa
from nmf_utilities import pivot_props
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker
from scipy import stats

from proteomics import get_pca_umap
from utils import plot_pca_ax, get_pcs

# --------------------------------------------------
# Color Schemes
# --------------------------------------------------
CPTAC_CMAP = {
    'GBM': '#52b788',
    'MEDUL': '#193e2e',
    'HNSCC': '#b7c47d',
    'LSCC': '#91bdff',
    'LUAD': '#1a759f',
    'CCRCC': '#ffd966',
    'COAD': '#ff8c69',
    'PDAC': '#962a13',
    'UCEC': '#5a189a',
    'OV': '#cdb4db',
    'BRCA': '#cd6090'
}

PAM50_SUBTYPE_CMAP = {
    "Basal":"#D81B60",
    "LumA": "#004D40",
    "LumB":"#1E88E5",
    "Her2":"#9250FE",
    "Normal-like":"#FFC107"
}

# OLD COLORS
# CPTAC_CMAP = {
#     'BRCA': (0.9677975592919913, 0.44127456009157356, 0.5358103155058701),
#     'LUAD': (0.5920891529639701, 0.6418467016378244, 0.1935069134991043),
#     'LSCC': (0.8087954113106306, 0.5634700050056693, 0.19502642696727285),
#     'MEDUL': (0.21044753832183283, 0.6773105080456748, 0.6433941168468681),
#     'UCEC': (0.6423044349219739, 0.5497680051256467, 0.9582651433656727),
#     'GBM': (0.33999999999999997, 0.8287999999999999, 0.86),
#     'COAD': (0.5529411764705883, 0.6274509803921569, 0.796078431372549),
#     'HNSCC': (0.7019607843137254, 0.7019607843137254, 0.7019607843137254),
#     'PDAC': (0.9882352941176471, 0.5529411764705883, 0.3843137254901961),
#     'CCRCC': (0.6509803921568628, 0.8470588235294118, 0.32941176470588235),
#     'OV': (0.03137254901960784, 0.28161476355247983, 0.5582622068435218)
# }

COHORT_CMAP = LinearSegmentedColormap('cohort_map', CPTAC_CMAP)

PTM_SCHEME = {
    'Protein': (0.9677975592919913, 0.44127456009157356, 0.5358103155058701),
    'Phosphorylation': (0.3126890019504329, 0.6928754610296064, 0.1923704830330379),
    'Acetylation': (0.23299120924703914, 0.639586552066035, 0.9260706093977744),
    #'RNA': (0.86, 0.3712, 0.33999999999999997)
    'RNA': (0.997341, 0.733545, 0.505167)
}

PTM_SCHEME2 = {
    'Proteome': (0.9677975592919913, 0.44127456009157356, 0.5358103155058701),
    'Phosphoproteome': (0.3126890019504329, 0.6928754610296064, 0.1923704830330379),
    'Acetylome': (0.23299120924703914, 0.639586552066035, 0.9260706093977744),
    'Transcriptome':(0.997341, 0.733545, 0.505167)
    #'Transcriptome': (0.86, 0.3712, 0.33999999999999997)
}

def _format_ax(ax):
    """Format axis."""
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    for axis in ['bottom','left',]:
        ax.spines[axis].set_linewidth(1.2)

# -------------------------
# QC
# -------------------------
def plot_cohort_coverage(full_feats_df, figsize=(18,8)):
    """
    Plot sample numbers & n-features detected.
    """
    fig,axes = plt.subplots(2,2,figsize=figsize)

    sns.barplot(
        x='Assay',
        y='samples',
        hue='Cohort',
        data=full_feats_df[full_feats_df['Type']=='Tumor'],
        palette=CPTAC_CMAP,
        edgecolor='black',
        ax=axes[0,0]
    )

    sns.barplot(
        x='Assay',
        y='samples',
        hue='Cohort',
        data=full_feats_df[full_feats_df['Type']!='Tumor'],
        palette=CPTAC_CMAP,
        edgecolor='black',
        ax=axes[0,1]
    )

    sns.barplot(
        x='Assay',
        y='coverage',
        hue='Cohort',
        data=full_feats_df[full_feats_df['Type']=='Tumor'],
        palette=CPTAC_CMAP,
        edgecolor='black',
        ax=axes[1,0]
    )

    sns.barplot(
        x='Assay',
        y='coverage',
        hue='Cohort',
        data=full_feats_df[full_feats_df['Type']!='Tumor'],
        palette=CPTAC_CMAP,
        edgecolor='black',
        ax=axes[1,1]
    )

    axes[0,0].set_ylabel("# Samples", fontsize=14)
    axes[0,1].set_ylabel("")
    axes[0,0].set_xlabel("")
    axes[0,1].set_xlabel("")

    axes[1,0].set_ylabel("Coverage", fontsize=14)
    axes[1,1].set_ylabel("")
    axes[1,0].set_xlabel("")
    axes[1,1].set_xlabel("")

    axes[0,0].set_xticklabels(['Proteome','Phospho','Acetyl'])
    axes[0,1].set_xticklabels(['Proteome','Phospho','Acetyl'])
    axes[1,0].set_xticklabels(['Proteome','Phospho','Acetyl'])
    axes[1,1].set_xticklabels(['Proteome','Phospho','Acetyl'])

    axes[0,0].set_title("Tumor", fontsize=16)
    axes[0,1].set_title("Normal", fontsize=16)

    axes[0,0].legend().remove()
    axes[1,0].legend().remove()
    axes[1,1].legend().remove()
    handles, labels = axes[0,1].get_legend_handles_labels()
    axes[0,1].legend(handles=handles, labels=labels, loc='center left', bbox_to_anchor=(1, .5), frameon=False)

    plt.tight_layout()

def plot_prot_qc(X, na_proportion_lim=0.7, title=None, axes=None, add_cbar=False, ylim=None):
    """
    Make masterplot.
    """
    X = X.copy().astype(float)

    proportion_lim = 100*na_proportion_lim

    # Aggr metrics
    X_agg = pd.DataFrame(100*(X.isna().sum(1) / X.shape[1]), columns=['p_missing'])
    X_agg['mu'] = X.median(1)
    X_agg['var'] = X.var(1)
    X_agg['keep'] = X_agg['p_missing']<proportion_lim

    # Create plot
    if axes is None:
        fig,ax = plt.subplots(figsize=(6,6))
        if add_cbar:
            cax = fig.add_axes([0.15, 0.2, 0.1, 0.02])
        ax2 = fig.add_axes([0.92,.125, 0.125, .755])
        ax3 = fig.add_axes([0.125,.9, 0.775, .125])
    else:
        ax,ax2,ax3 = axes
        if add_cbar:
            cax = fig.add_axes([0.15, 0.2, 0.1, 0.02])

    ax.scatter(
        X_agg[~X_agg['keep']]['p_missing'],
        X_agg[~X_agg['keep']]['mu'].fillna(0),
        c='lightgrey',
        alpha=0.1,
        vmax=10,
        rasterized=True
    )
    im = ax.scatter(
        X_agg[X_agg['keep']]['p_missing'],
        X_agg[X_agg['keep']]['mu'],
        c=X_agg[X_agg['keep']]['var'],
        alpha=0.5,
        vmax=10,
        vmin=0,
        cmap='coolwarm',
        rasterized=True
    )

    ax.axhline(0, c='black', linestyle=':')

    if add_cbar:
        cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
        cbar.set_ticks([0,5,10])
        cax.set_title(r'$\sigma^2$')

    if ylim is None:
        ax.set_ylim(-np.abs(X_agg['mu']).max()*1.05,np.abs(X_agg['mu']).max()*1.05)
    else:
        ax.set_ylim(ylim)

    ax.set_xlim(-5,105)

    # Mean Dist
    sns.kdeplot(X_agg[~X_agg['keep']]['mu'], vertical=True, ax=ax2, color='grey', label=None)
    sns.kdeplot(X_agg[X_agg['keep']]['mu'], vertical=True, ax=ax2, color='darkblue', label=None)
    ax2.set_yticks([])
    ax2.set_xticks([])
    ax2.set_ylim(*ax.get_ylim())
    ax2.axhline(0, c='black', linestyle=':')
    ax2.legend().remove()

    # P dist
    sns.kdeplot(X_agg['p_missing'], ax=ax3, color='darkblue', shade=True, label=None, cut=0)
    ax3.set_yticks([])
    ax3.set_xticks([])
    ax3.set_xlim(*ax.get_xlim())
    ax3.axvline(proportion_lim, c='black', linestyle=':')
    ax3.legend().remove()

    ax.set_xlabel("Missing %", fontsize=14)
    ax.set_ylabel("Median-MAD LogFC", fontsize=14)
    ax.set_xticks([0,25,50,75,100])

    n_sites_remove = X_agg[~X_agg['keep']].shape[0]
    n_sites_keep = X_agg[X_agg['keep']].shape[0]

    ax.text(0,ax.get_ylim()[1]*.8,'$m = {}$'.format(n_sites_keep), ha='left', fontsize=12)
    ax.text(100,ax.get_ylim()[1]*.8,'$m = {}$'.format(n_sites_remove), ha='right', fontsize=12)

    ax.axvline(proportion_lim, c='black', linestyle=':')
    plt.title(title, fontsize=16)

    return X_agg

def plot_prot_grid_qc(raw_data_files, cohorts, assays, lim):
    """
    Create grid QC plot.
    """
    import matplotlib.gridspec as gridspec

    fig = plt.figure(figsize=(len(assays)*8,cohorts.shape[0]*8))
    gs = gridspec.GridSpec(9*cohorts.shape[0], 10*len(assays), hspace=0.5, wspace=0.5)

    aggr_df = list()

    for idx,cohort in enumerate(cohorts):
        for col_idx,assay in enumerate(assays):
            ax = plt.subplot(gs[(idx)*8+2:(idx+1)*8-1, 10*col_idx:(10*col_idx+7)])
            ax2 = plt.subplot(gs[(idx)*8+2:(idx+1)*8-1,(10*col_idx+7):(10*col_idx+9)])
            ax3 = plt.subplot(gs[(idx)*8:(idx)*8+2, 10*col_idx:(10*col_idx+7)])

            try:
                _aggr_df = plot_prot_qc(
                    raw_data_files[assay][cohort]['Tumor'],
                    na_proportion_lim=lim[col_idx],
                    title='{} {}'.format(cohort, assay.capitalize()),
                    axes=(ax,ax2,ax3),
                    ylim=(-10,10)
                )

                ax.set_yticks([-10,-5,0,5,10])

                _aggr_df['cohort'] = cohort
                _aggr_df['assay'] = assay
                aggr_df.append(_aggr_df)

            except:
                print(assay, cohort)

    return aggr_df

def plot_prot_upset(features_df, assay, tn, pmap_df=None, min_subset_size=50, order=None, ncols=20):
    """
    Create upset plot.
    """
    import upsetplot as upsetplot
    from upsetplot import UpSet

    bc_counts_df = features_df[features_df['type']==tn].groupby([assay,'cohort']).size().reset_index().pivot(index=assay, columns='cohort', values=0)
    bc_counts_df = bc_counts_df > 0

    if order is not None:
        bc_counts_df = bc_counts_df[order]
        sort_categories_by = None
    else:
        sort_categories_by="cardinality"

    # Cohorts
    cohorts = bc_counts_df.columns

    # Add PTM-Sites per protein
    if assay == 'proteome':
        assert pmap_df is not None, "Provide protein map."
        phosph_s = pmap_df[pmap_df['feature']=='phosphoproteome'].groupby('accession_number').size()
        phosph_s.name = "# Phospho"

        acetyl_s = pmap_df[pmap_df['feature']=='acetylome'].groupby('accession_number').size()
        acetyl_s.name = "# Acetyl"

        bc_counts_df = bc_counts_df.join(phosph_s).join(acetyl_s).fillna(0)

    bc_counts_df = bc_counts_df.join(pmap_df[['protein_mw','accession_number']]).rename(columns={'protein_mw':'MW'})

    # Index by cohorts
    if ncols is not None:
        gb = bc_counts_df.groupby(list(cohorts)).size().sort_values(ascending=False)
        min_subset_size = gb.iloc[ncols-1]

    bc_full_df = bc_counts_df.copy()
    bc_counts_df = bc_counts_df.set_index(list(cohorts))

    upset = UpSet(
        bc_counts_df,
        sort_by='cardinality',
        min_subset_size=min_subset_size,
        sort_categories_by = sort_categories_by,
        show_counts=False,
        totals_cdict=CPTAC_CMAP
    )

    for c in cohorts:
        upset.style_subsets(present=c, edgecolor="black", facecolor='darkgrey', linewidth=1)
        #upset.style_subsets(present=c, absent=[x for x in cohorts if x!=c], edgecolor="black", facecolor=CPTAC_CMAP[c], linewidth=1)

    if assay == 'proteome':
        #upset.add_catplot(value="# Phospho", kind='box', color='darkgrey')
        #upset.add_catplot(value="# Acetyl", kind='box', color='darkgrey')
        upset.add_catplot(value="MW", kind='point', color='darkgrey',  ylog_scale=True,)
    elif assay == 'phosphoproteome' or assay=='acetylome':
        upset.add_catplot(value="MW", kind='point', color='darkgrey', drop_index_duplicate=True, ylog_scale=True, s=2, alpha=0.3)
        #upset.add_catplot(value="MW", kind='box', color='darkgrey', drop_index_duplicate=True, ylog_scale=True)

    upset.plot()
    return bc_full_df

def plot_pca_grid(X, cohort_s=None, cohort_colors=None, normalize=True, C=None, title=None, retFig=False):
    """
    plot PCA grid.
    """
    if normalize is False:
        X = X - X.mean(axis=0)
        X = X / np.sqrt(X.pow(2).sum(axis=0))

    P_df, pca, gns = get_pcs(X, normalize=normalize, return_genes=True, C=C)
    fig,axes = plt.subplots(2,2,figsize=(10,10))

    cohort_s = cohort_s.loc[P_df.index]

    plot_pca_ax(P_df, pca, ax=axes[0,0], cohort_s=cohort_s, cohort_colors=cohort_colors,order=[1,2,3])
    plot_pca_ax(P_df, pca, ax=axes[0,1], cohort_s=cohort_s, cohort_colors=cohort_colors, order=[1,4,4])
    plot_pca_ax(P_df, pca, ax=axes[1,0], cohort_s=cohort_s, cohort_colors=cohort_colors, order=[3,2,5])
    plot_pca_ax(P_df, pca, ax=axes[1,1], cohort_s=cohort_s, cohort_colors=cohort_colors, order=[3,4,6])

    axes[0,0].legend().remove()
    axes[0,1].legend(loc='center right', bbox_to_anchor=(1.35, 0))
    axes[1,0].legend().remove()
    axes[1,1].legend().remove()

    if title is not None:
        fig.suptitle(title,y=0.925, fontsize=18)

    if retFig:
        return fig, axes

# --------------------------------------------------
# Plotting Functions
# --------------------------------------------------

def plot_proportions(df: pd.DataFrame, figsize: tuple = (6,8), color=None, ax=None):
    """
    Plot Proportions
    ------------------
    Args:
        * df: pd.DataFrame of (clusters x cohorts)
        * figsize
        * cmap: matplotlib colormap to use
        * ax: axis

    Returns:
        * fig
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    df.iloc[::-1].plot(
        kind='barh',
        ax=ax,
        stacked=True,
        width=1,
        edgecolor='black',
        color=color,
        rasterized=True
    )

    if np.max(np.max(df)) <= 1.0:
        ax.set_xlim(0,1)

    ax.set_ylim(ax.get_ylim()[0]+0.25, ax.get_ylim()[1]-0.25)
    ax.set_xlabel(df.columns.name.capitalize(), fontsize=20)
    ax.legend(bbox_to_anchor=(1, 1.0125))
    ax.set_ylabel(df.index.name.capitalize(), fontsize=20)

def plot_site_dist(df, attr, ax=None, diff=0.5, max_norm=0.5):
    """
    Plot site distribution from signatures.
    --------------------------
    Args:
        * df: pd.DataFrame (sites x signature information)
        * attr: variable to use
        * ax: axis
        * diff: filter for sigantures
        * max_norm: filter for signatures
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(4,4))

    df = df[(df['diff']>diff) & (df['max_norm']>max_norm)]

    for feat in np.unique(df['index']):
        sns.distplot(
            df[df['index']==feat][attr],
            ax=ax,
            label="{} (n={})".format(feat, df[df['index']==feat].shape[0]),
            hist=False,
            kde_kws={"linewidth":2.5},
            color=PTM_SCHEME[feat]
        )

    ax.legend(frameon=False)
    ax.set_xlim(np.min(df[attr]))

def plot_ptm_prot_ols(ptm_df, prot_df, meta, cohort, aggr=None, ax=None, l=None, plot_density=False, collapsed=False, ptm2Prot_dict=None):
    """
    Plot PTM Protein OLS
    -----------------------
    Args:
        * ptm: ptm data
        * prot: proteomics data
        * meta: metadata
        * cohort: cohort to use
        * aggr: how to aggregate (default None)
        * ax
        * l axis square limit

    Return:
        * pd.DataFrame with residuals
        * dict: mapping cohort -> statsmodels fit
    """
    import statsmodels.api as sm
    from tqdm import tqdm
    import proteomics as prot

    if ax is None:
        fig,ax = plt.subplots(figsize=(5,5))

    regression_dict = {}
    df = list()

    ptm_df.index.name = 'id'
    prot_df.index.name = 'id'

    df = prot.ptm_comparison(ptm_df[meta[meta==cohort].index], prot_df[meta[meta==cohort].index], aggr=aggr, collapsed=collapsed, ptm2Prot_dict=ptm2Prot_dict)
    df['cohort'] = cohort
    df = df.dropna(subset=('PTM','Proteome'))

    if l is None:
        amax_ptm = np.max(np.abs(np.array([np.max(df['PTM']),np.min(df['PTM'])])))
        amax_prot = np.max(np.abs(np.array([np.max(df['Proteome']),np.min(df['Proteome'])])))
        l = np.max([amax_ptm,amax_prot])

    ax.set_xlim(-l,l)
    ax.set_ylim(-l,l)

    ax.axhline(0, color='k', alpha=0.1, zorder=-1)
    ax.axvline(0, color='k', alpha=0.1, zorder=-1)

    if plot_density:
        hexdensity(
            df['Proteome'],
            df['PTM'],
            cmap='viridis',
            scale=None,
            unit='Med-MAD LFC',
            ax=ax,
            entity='sites',
            rasterized=True
        )
    else:
        sns.regplot(
            x="PTM",
            y="Proteome",
            data=df,
            color=CPTAC_CMAP[cohort],
            ax=ax,
            ci=68,
            scatter_kws={'alpha':0.05, 'rasterized':True},
            truncate=False
        )

    print("   * {} / {} sites with matching protein".format(np.unique(df.index).shape[0], ptm_df.index.shape[0]))

def plot_signature_2x2_grid(H, meta_s, signatures, site_annot_s, diff=0.5, max_norm=0.5, figsize=(9,8)):
    """
    Plot Signature 2x2 Grid
    -------------------------
    Args:
        * H: H matrix from NMF
        * meta_s: cohort series
        * signatures: signatures from signatureanalyzer
        * site_annot_s: site annotation series of PTM/Protein type
        * diff: signature filter
        * max_norm: signature filter
        * figsize
    """
    fig, axes = plt.subplots(2,2,figsize=figsize)

    _color = [CPTAC_CMAP[x] for x in (np.unique(meta_s))]
    plot_proportions(pivot_props(H[['max_id']].join(meta_s), norm=False, y='max_id'), color=_color, ax=axes[0,0])
    plot_proportions(pivot_props(H[['max_id']].join(meta_s), norm=True, y='max_id'), color=_color, ax=axes[0,1])
    plot_site_dist(signatures.join(site_annot_s), 'diff', ax=axes[1,0], diff=diff, max_norm=max_norm)
    plot_site_dist(signatures.join(site_annot_s), 'max_norm', ax=axes[1,1], diff=diff, max_norm=max_norm)

    axes[0,0].legend().remove()
    axes[0,0].set_xlabel('')
    axes[0,0].set_ylabel('Signature', fontsize=12)
    axes[0,1].set_xlabel('')
    axes[0,1].set_ylabel('')
    axes[1,0].legend().remove()

    axes[1,0].set_xlabel("Norm. Factor Difference", fontsize=12)
    axes[1,1].set_xlabel("Norm. Factor Weight", fontsize=12)

    plt.tight_layout()
    return fig

def plot_signature_2x1_grid(df, meta_s, cluster_idx='clusters', cluster_name='Consensus Cluster', figsize=(9,4)):
    """
    Plot Signature 2x1 Grid
    -------------------------
    Args:
        * H: H matrix from NMF
        * meta_s: cohort series
        * cluster_idx: index of input dataframe to cluster by
        * cluster_name: y-axis name of plot
        * figsize: tuple (int,int)

    Returns:
        * fig
    """
    fig, axes = plt.subplots(1,2,figsize=figsize)

    _color = [CPTAC_CMAP[x] for x in np.unique(list(CPTAC_CMAP.keys()))]
    plot_proportions(pivot_props(df[[cluster_idx]].join(meta_s), norm=False, y=cluster_idx), color=_color, ax=axes[0])
    plot_proportions(pivot_props(df[[cluster_idx]].join(meta_s), norm=True, y=cluster_idx), color=_color, ax=axes[1])

    axes[0].legend().remove()
    axes[0].set_xlabel('Counts', fontsize=12)
    axes[0].set_ylabel(cluster_name, fontsize=12)
    axes[1].set_xlabel('Proportion', fontsize=12)
    axes[1].set_ylabel('')

    plt.tight_layout()
    return fig

def plot_proteomics_heatmap(
    X: pd.DataFrame,
    signatures: pd.DataFrame,
    order_series: pd.Series,
    ptm: Union[pd.Series,None] = None,
    diff: float = 0.5,
    max_norm: float = 0.5,
    figsize: tuple = (16,12),
    cmap: str ="YlGnBu",
    display_y: bool = False,
    vmax: float = None,
    vmin: float = None,
    cohort_s: Union[pd.Series,None] = None,
    y_hm_label: str = 'Genes',
    cbar_hm_label: str = 'Normalized Expression',
    metas: Union[None,list] = None
    ):
    """
    Plot marker map.
    -----------------------------
    Args:
        * X: pd.DataFrame of input sample x feature matrix
        * signatures: pd.DataFrame signatures output;
            this bundles information about the weightings of each feature (ex. gene) and
            what signature they map to
        * order_series: series of samples mapping to subgroups
            index: X.index
            values: subgrouping
        * subset_series: a pd.Series with the index as the gene name or ID that
            matches the marker matrix & has a "Subgroup" column for labeling
        * diff: difference of loading for called signature vs. rest
        * max_norm: strength of loading for called signature
        * figsize: size of figure
        * cmap: colormap for plot
        * display_y: whether or not to display feature names
        * vmax: colorbar max
        * vmin: colorbar min
        * cohort_s: cohort_series dataframe (added on top of plot)
        * y_hm_label: label of y-axis on heatmap (ex. Genes, Protein LFC, etc.)
        * cbar_hm_label: label of heatmap colorbar

    Returns:
        * plt.Figure
    """
    # Filter for marker PTMS
    signatures_filt = signatures[(signatures['diff'] > diff) & (signatures['max_norm'] > max_norm)]

    # Remove signatures with no marker genes associated
    order_series = order_series[order_series.isin(set(signatures_filt['max_id'].astype(int)))]

    sig_index_list = []
    for s in np.unique(order_series.sort_values().values):
        x = signatures_filt[signatures_filt['max_id']==s]
        x['pos'] = x.reset_index()['index'].apply(lambda x: x[-1:]!='n').values
        sig_index_list += list(x.sort_values(['pos','diff',], ascending=False).index)

    # Filter X matrix
    sample_markers = X.loc[sig_index_list, order_series.sort_values().index]

    # Set horizontal lines
    hz_lines = np.unique(sample_markers.join(signatures_filt).loc[:,'max_id'].values, return_index=True)[1]

    fig, ax = plt.subplots(figsize=figsize)

    x0 = ax.get_position().x0
    x1 = ax.get_position().x1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1
    buf = y1*0.01

    cbar_ax = fig.add_axes([.91, 0.5, .025, .3])

    sns.heatmap(sample_markers, ax=ax, cmap=cmap, rasterized=True, cbar_ax=cbar_ax, vmax=vmax, vmin=vmin)
    v,c = np.unique(order_series, return_counts=True)

    # plot horizontal lines
    _c = np.cumsum(c)
    _ci = np.roll(_c,2)
    _ci[0] = 0
    _ci[1] = 0
    ax.hlines(hz_lines, _ci, _c, rasterized=True)

    # plot vertical lines
    _h = list(hz_lines)
    _h.append(sample_markers.shape[0])
    ax.vlines(np.cumsum(c)[:-1], _h[:-2], _h[2:], rasterized=True)
    ax.vlines(np.cumsum(c)[:-1], 0, sample_markers.shape[0], alpha=0.4, rasterized=True)

    # set ticks
    ax.set_xticks(np.cumsum(c)-c/2)
    ax.set_xticklabels(v, rotation=360,fontsize=14)
    ax.set_yticks(np.arange(sample_markers.index.values.shape[0]))

    # add gene markings
    if ptm is not None:
        ax.set_yticks([])
        ax.set_yticklabels([], rasterized=True)

        cdict = PTM_SCHEME

        if "RNA" in ptm.values:
            order_dict = {PTM_SCHEME['Protein']:0, PTM_SCHEME['Phosphorylation']:1, PTM_SCHEME['RNA']: 3}
        else:
            order_dict = {PTM_SCHEME['Protein']:0, PTM_SCHEME['Phosphorylation']:1, PTM_SCHEME['Acetylation']:2}

        lax = fig.add_axes([x0-4*buf, y0, 2*buf, y1-y0])
        lax.set_xticks([])
        lax.set_yticks([])

        meta = sample_markers.drop(columns=list(sample_markers)).join(ptm).iloc[:,0]

        colors_conversion, meta_colormap = sa.pl.series_to_colors(meta, cdict=cdict)
        meta_colormap_inv = dict([[v,k] for k,v in meta_colormap.items()])
        meta_colormap_inv = {(k[0],k[1],k[2]):v for k,v in meta_colormap_inv.items()}
        cbar_lax = fig.add_axes([x0-8*buf, y1-x1*.25*.75, 2*buf, x1*.25*.75])

        # Add heatmapping
        mat,cmap = sa.pl.color_list_to_matrix_and_cmap(colors_conversion, order_dict=order_dict)
        sns.heatmap(
            mat.T,
            cmap=cmap,
            ax=lax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=cbar_lax,
        )

        cb_ticks = [float(t.get_text().replace('−','-')) for t in cbar_lax.get_yticklabels()]

        color_value_mapping = dict()
        for v in np.unique(mat):
            color_code = list(cmap.__call__(v))
            color_code = tuple(color_code[:3])
            color_value_mapping[v] = meta_colormap_inv[color_code]

        cbar_lax.get_yaxis().set_ticks([])

        n_labels = len(list(color_value_mapping.keys()))
        vals = [x * ((n_labels)/(n_labels+1)) + 0.5 * ((n_labels)/(n_labels+1)) for x in list(color_value_mapping.keys())]

        if "RNA" in ptm.values:
            #print(vals)
            #cbar_lax.get_yaxis().set_ticks([0])
            pass
        else:
            cbar_lax.get_yaxis().set_ticks([0.35, 1, 1.65])

        cbar_lax.get_yaxis().set_ticklabels(list(color_value_mapping.values()),)
        cbar_lax.yaxis.set_ticks_position('left')

        cbar_lax.set_frame_on(True)

        lax.set_ylabel("")
        ax.set_ylabel("")

        for _, spine in lax.spines.items():
            spine.set_visible(True)

    # --------------subset genes-------------------
    else:
        if display_y:
            ax.set_yticklabels(sample_markers.index.values, fontsize=5, rasterized=True)
        else:
            ax.set_yticks([])
            ax.set_yticklabels([], rasterized=True)

        ax.set_ylabel(y_hm_label, fontsize=14)

    upper_buffer = y1+2*buf

    # Add sample annot
    if metas is not None:
        for idx,meta in enumerate(metas):
            new_ax = [x0, y1+(idx+2)*buf+(idx*2)*buf, x1*.861, 2*buf]
            lax = fig.add_axes(new_ax)

            if meta.name == 'cohort':
                cdict = CPTAC_CMAP
                order_dict = None
                #order_dict = {COLOR_SCHEME['pink']:0, COLOR_SCHEME['yellow']:1, COLOR_SCHEME['teal']:2, COLOR_SCHEME['purple']:3,  COLOR_SCHEME['orange']:4,  COLOR_SCHEME['lightblue']:5}
                #order_dict = {COLOR_SCHEME['pink']:0, COLOR_SCHEME['yellow']:1, COLOR_SCHEME['teal']:2, COLOR_SCHEME['purple']:3}
            elif meta.name == 'grading':
                cdict = TUMOR_STAGE_SCHEME
                order_dict = None
            else:
                cdict = None
                order_dict = None

            meta = meta.loc[sample_markers.columns]
            _meta_v = pd.unique(meta)

            if None in _meta_v and True in _meta_v and len(_meta_v)==2:
                cdict = {True:(0.21044753832183283, 0.6773105080456748,0.6433941168468681), None:(1.0,1.0,1.0,1.0)}
                order_dict = {(1.0,1.0,1.0,1.0):0,(0.21044753832183283,0.6773105080456748,0.6433941168468681):1}

            cluster_color_list, _ = sa.pl.series_to_colors(meta,cdict=cdict)
            mat,cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list,order_dict=order_dict)
            sns.heatmap(mat, cmap=cmap, ax=lax, yticklabels=False, xticklabels=False, cbar=False)
            lax.set_ylabel(meta.name, y=0, ha='right', rotation=360)

            for _, spine in lax.spines.items():
                spine.set_visible(True)

        upper_buffer = y1+(idx+3)*buf+(idx*2)*buf

    if cohort_s is not None:
        # Get ordering and samples
        cohort_s = cohort_s.loc[sample_markers.columns]
        cdict = CPTAC_CMAP
        order_dict = None

        # Create axis
        cs_ax = fig.add_axes([x0, y1+2*buf, x1*.861, 2*buf])
        cs_ax.set_xticks([])
        cs_ax.set_yticks([])

        cbar_cs_ax = fig.add_axes([x0, y1+6*buf, x1*.25, 2*buf])

        colors_conversion, meta_colormap = sa.pl.series_to_colors(cohort_s, cdict=cdict)
        meta_colormap_inv = dict([[v,k] for k,v in meta_colormap.items()])
        meta_colormap_inv = {(k[0],k[1],k[2]):v for k,v in meta_colormap_inv.items()}
        mat,cmap = sa.pl.color_list_to_matrix_and_cmap(colors_conversion, order_dict=order_dict)

        sns.heatmap(
            mat,
            cmap=cmap,
            ax=cs_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=cbar_cs_ax,
            cbar_kws={"orientation": "horizontal"}
        )

        cb_ticks = [float(t.get_text().replace('−','-')) for t in cbar_cs_ax.get_yticklabels()]
        cbar_cs_ax.xaxis.set_ticks([])

        cbar_cs_ax.set_frame_on(True)
        [spine.set_visible(True) for _, spine in cs_ax.spines.items()]
    # --------------sample annot-------------------

    ax.set_title('')
    ax.set_xlabel('NMF Signatures', fontsize=14)

    cbar_ax.set_ylabel(cbar_hm_label, fontsize=12)

    [spine.set_visible(True) for _, spine in ax.spines.items()]

    return fig

def plot_sig_props(H, meta_s, signatures, figsize=(9,8)):
    """
    Plot Sample Signature Proportions by Meta
    -------------------------
    Args:
        * H: H matrix from NMF
        * meta_s: cohort series
        * signatures: signatures from signatureanalyzer
        * site_annot_s: site annotation series of PTM/Protein type
        * diff: signature filter
        * max_norm: signature filter
        * figsize
    """
    import matplotlib.gridspec as gridspec

    fig = plt.figure(constrained_layout=True, figsize=figsize)
    spec = gridspec.GridSpec(ncols=5, nrows=5, figure=fig)

    ax1 = fig.add_subplot(spec[:, :4])
    ax2 = fig.add_subplot(spec[:, 4:])

    _color = [CPTAC_CMAP[x] for x in (np.unique(meta_s))]
    plot_proportions(pivot_props(H[['max_id']].join(meta_s), norm=False, y='max_id'), color=_color, ax=ax1)
    plot_proportions(pivot_props(H[['max_id']].join(meta_s), norm=True, y='max_id'), color=_color, ax=ax2)

    ax1.set_xlabel("# of Samples", fontsize=24)
    ax1.legend().remove()
    ax1.set_ylabel("")

    ax2.set_xlabel("Proportion", fontsize=24)
    ax2.set_ylabel("")

    ax1.set_yticklabels(ax1.get_yticklabels(), fontsize=18)
    ax2.set_yticks([])

    ax2.legend(fontsize=16, loc=(1.2,.66), frameon=False)

    plt.tight_layout()
    return fig

def plot_sig_W_props(W, pmap_df, figsize=(10,12)):
    """
    Plot Feature Signature Proportions by Assay
    -------------------------
    Args:
        * W: W-matrix from NMF
        * pmap_df: mapping of features to assay types
    """
    df = W[['max_id','max','max_norm']]
    df['n'] = [-1 if x else 1 for x in df.index.str.endswith('_n')]
    df['max'] = df['max']*df['n']
    df['max_norm'] = df['max_norm']*df['n']
    df.index = df.index.str.replace('_n','')
    df = df.join(pmap_df[['geneSymbol','feature']])
    df['qval'] = 0
    df['feature'] = df['feature'].str.capitalize()

    fig,axes = plt.subplots(1, 2, figsize=figsize, sharey=True)

    plot_proportions(pivot_props(df[df['n']==-1], norm=False, x='feature',y='max_id'), color=PTM_SCHEME2, ax=axes[0])
    plot_proportions(pivot_props(df[df['n']==1], norm=False, x='feature',y='max_id'), color=PTM_SCHEME2, ax=axes[1])

    axes[0].legend().remove()
    axes[1].legend().remove()
    axes[0].set_xlim(1100,0)
    axes[1].set_xlim(0,1100)
    axes[0].set_title("Down-Regulated", fontsize=20, x=1, ha='right')
    axes[1].set_title("Up-Regulated", fontsize=20, x=0, ha='left')
    axes[0].yaxis.tick_right()
    axes[0].set_yticklabels(axes[0].get_yticklabels(), fontsize=16, ha='left')
    axes[0].set_xlabel("Features per Signature", fontsize=18)
    axes[1].set_xlabel("Features per Signature", fontsize=18)
    axes[0].set_ylabel("")
    axes[0].spines['left'].set_visible(False)
    axes[0].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].spines['top'].set_visible(False)

    plt.tight_layout()

def plot_qc_rna_grid(
    metrics,
    c_s,
    add_legend=True,
    figsize=(8,8),
    format_ax=True,
    series_cdict=None,
    **legend_kwargs
    ):
    """
    Plot QC grid.
    """
    import signatureanalyzer as sa
    color_s, color_dict = sa.pl.series_to_colors(c_s, cdict=series_cdict)

    # Plot figure
    fig,axes = plt.subplots(2,2,figsize=(8,8), sharey=True)

    def _format_ax(ax):
        if format_ax:
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            for axis in ['bottom','left',]:
                ax.spines[axis].set_linewidth(1.2)

    # --------- 0,0 ---------
    cax = fig.add_axes([0.4, 0.925, 0.1, 0.01])
    im = axes[0,0].scatter(
        metrics["Duplicate Rate of Mapped"],
        metrics['Genes Detected'],
        s=20,
        alpha=0.8,
        c=np.log(metrics['Unique Mapping, Vendor QC Passed Reads']),
        rasterized=True
    )

    fig.colorbar(im, cax=cax, orientation='horizontal')
    cax.set_title("log Unique Mapping", fontsize=8)

    axes[0,0].set_ylabel('Genes Detected',fontsize=16)
    _format_ax(axes[0,0])

    # --------- 0,1 ---------
    im = axes[0,1].scatter(
        metrics["Median Exon CV"],
        metrics['Genes Detected'],
        s=20,
        alpha=0.8,
        c=np.log(metrics['Unique Mapping, Vendor QC Passed Reads']),
        rasterized=True
    )
    axes[0,1].set_title("n={}".format(metrics.shape[0]), fontsize=16, x=0.85, y=0.85, fontstyle='italic')
    _format_ax(axes[0,1])

    # --------- 1,0 ---------
    axes[1,0].scatter(
        metrics["Duplicate Rate of Mapped"],
        metrics['Genes Detected'],
        s=20,
        alpha=0.8,
        c=color_s,
        rasterized=True
    )

    axes[1,0].set_ylabel('Genes Detected',fontsize=16)
    axes[1,0].set_xlabel('Duplicate Rate',fontsize=16)
    _format_ax(axes[1,0])

    # --------- 1,1 ---------
    axes[1,1].scatter(
        metrics["Median Exon CV"],
        metrics['Genes Detected'],
        s=20,
        alpha=0.8,
        c=color_s,
        rasterized=True
    )
    axes[1,1].set_xlabel('Median Exon CV',fontsize=16)
    _format_ax(axes[1,1])

    axes[1,1].set_xlim(*axes[1,1].get_xlim())
    axes[1,1].set_ylim(*axes[1,1].get_ylim())

    plt.tight_layout()

    if add_legend:
        for k,v in color_dict.items():
            m = axes[1,1].scatter(-1e2, -1e2, alpha=.8, c=v, label=k)

        axes[1,1].legend(**legend_kwargs)

def plot_n_top_hvg(tpm_norm, gene_name, n_top=2500, mean_thresh=None):
    """
    Plot n top hvg.
    """
    if mean_thresh is not None:
        x = tpm_norm.mean(1)
        prev_g = tpm_norm.shape[0]
        tpm_norm = tpm_norm.loc[x[x>mean_thresh].index]
        print("filtered out {} / {} lowly expressed genes.".format(prev_g-tpm_norm.shape[0], prev_g))

    gene_cv = tpm_norm.std(1) / tpm_norm.mean(1)
    n_top_genes = n_top

    fig,axes = plt.subplots(2,1, figsize=(10,6), sharex=True)

    lim = gene_cv.sort_values(ascending=False).iloc[n_top_genes]
    top_genes = gene_cv.sort_values(ascending=False).iloc[:n_top_genes]

    # coloring
    colors = {True:'red', False:'lightgrey'}
    _df = pd.DataFrame(gene_cv).join(gene_name)
    _df['keep'] = _df[0].index.isin(top_genes.index)

    sns.distplot(gene_cv, ax=axes[0], bins=100, kde=False, hist_kws=dict(edgecolor="k", linewidth=1,))
    axes[0].axvspan(0, _df[_df['keep']].sort_values(0)[0][0], zorder=0, alpha=0.2, color='lightgrey')
    axes[0].set_title('Selecting {} / {} Highly Variable Genes'.format(n_top_genes, gene_cv.shape[0]), fontsize=16)
    axes[0].set_ylabel("Frequency", fontsize=14)

    # Scatter plot
    axes[1].scatter(gene_cv, tpm_norm.mean(1),s=3, alpha=0.5, c=_df['keep'].apply(lambda x: colors[x]), rasterized=True)
    axes[1].set_xlabel("CV Across Samples", fontsize=14)
    axes[1].set_ylabel("Mean expression", fontsize=14)
    axes[1].set_xlim(0, axes[1].get_xlim()[1])
    axes[1].set_ylim(0, axes[1].get_ylim()[1])

    plt.tight_layout()

    return top_genes

def plot_sample_dist(X, meta_s, figsize=(10,4), ax=None, title=None, xlim=(-10,10), color_by='cohort'):
    """
    Plot sample distribution - for each sample, it plots the distribution
    of all sites. After median normalization and scaling, this should be
    zero-mean gaussian.
    -------------------------
    Args:
        * X: dict with each feature as a key and matrices
        * meta_s: pd.series of cohort
        * figsize: size of figure
        * ax: matplotlib axis
        * title: title for plot
        * xlim: x limits

    Returns:
        None
    """
    random_idx = X.sample(frac=1).index
    meta_s = meta_s.loc[random_idx].copy()
    X = X.loc[random_idx].copy()

    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    if color_by=='cohort':
        for idx in range(X.shape[0]):
            sns.distplot(
                X.iloc[idx],
                hist=False,
                color=CPTAC_CMAP[meta_s[idx]],
                ax=ax,
                kde_kws={'linewidth':1},
                label=None
            )
    else:
        sample_colors, color_map_dict = sa.pl.series_to_colors(meta_s)
        for idx in range(X.shape[0]):
            sns.distplot(
                X.iloc[idx],
                hist=False,
                ax=ax,
                kde_kws={'linewidth':1},
                color=color_map_dict[meta_s[idx]]
            )

    ax.set_xlim(xlim)
    ax.set_xlim(*ax.get_xlim())
    ax.set_ylim(*ax.get_ylim())

    if color_by=='cohort':
        for k in CPTAC_CMAP.keys():
            if k in pd.unique(meta_s):
                m = ax.scatter(-1, -1, alpha=.8, c=CPTAC_CMAP[k], label=k, rasterized=True)
    else:
        for k in color_map_dict.keys():
            if k in pd.unique(meta_s):
                m = ax.scatter(-1, -1, alpha=.8, c=np.array(color_map_dict[k])[np.newaxis,:], label=k, rasterized=True)

    ax.legend(frameon=False)
    ax.set_xlabel("")
    ax.axvline(0, color='black', linestyle=':', alpha=0.4)

    if title is not None:
        ax.set_title(title, fontsize=16, x=0.01, y=0.875, ha='left', weight='bold')

def plot_signature_dist(X_s, meta_s, figsize=(10,4), ax=None, title=None, ylim=None, color_by='cohort'):
    """
    Plot sample distribution - for each sample, it plots the distribution
    of all sites. After median normalization and scaling, this should be
    zero-mean gaussian.
    -------------------------
    Args:
        * X: dict with each feature as a key and matrices
        * meta_s: pd.series of cohort
        * figsize: size of figure
        * ax: matplotlib axis
        * title: title for plot
        * xlim: x limits

    Returns:
        None
    """
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    meta_s = meta_s.loc[X_s.index]
    if color_by=='cohort':
        for cohort in np.unique(meta_s):
            sns.distplot(
                X_s.loc[meta_s[meta_s==cohort].index],
                hist=False,
                color=CPTAC_CMAP[cohort],
                ax=ax,
                kde_kws={'linewidth':1},
                label=None
            )
    else:
        sample_colors, color_map_dict = sa.pl.series_to_colors(meta_s)
        for idx in range(X.shape[0]):
            sns.distplot(
                X.iloc[idx],
                hist=False,
                ax=ax,
                kde_kws={'linewidth':1},
                color=color_map_dict[meta_s[idx]]
            )

    ax.set_xlim(0, ax.get_xlim()[1])
    if ylim is not None:
        ax.set_ylim(*ylim)

    ax.set_ylim(*ax.get_ylim())

    if color_by=='cohort':
        for k in CPTAC_CMAP.keys():
            if k in pd.unique(meta_s):
                m = ax.scatter(-1, -1, alpha=.8, c=CPTAC_CMAP[k], label=k, rasterized=True)
    else:
        for k in color_map_dict.keys():
            if k in pd.unique(meta_s):
                m = ax.scatter(-1, -1, alpha=.8, c=np.array(color_map_dict[k])[np.newaxis,:], label=k, rasterized=True)

    ax.legend(frameon=False)
    ax.set_xlabel("")

    if title is not None:
        ax.set_title(title, fontsize=16, x=0.01, y=0.875, ha='left', weight='bold')

def plot_cm(
    s1,
    s2,
    figsize=(5,4),
    normalize=False,
    use_unique=False,
    cmap="Blues",
    linewidth=0.25,
    linecolor='lightgrey',
    **kwargs):
    """
    Plot simple confusion matrix.

    Args:
        * s1: series1
        * s2: series2
        * figsize: figure size
        * normalize: normalize by s2
            Normalizes dataframe
    """
    from sklearn.metrics import confusion_matrix

    labels = np.union1d(s1.astype('category').cat.categories, s2.astype('category').cat.categories)
    cm = pd.DataFrame(data=confusion_matrix(s1,s2), index=labels, columns=labels)

    if normalize:
        cm = cm / cm.sum(1)[:,np.newaxis]
        missing_left=True

    fig,ax = plt.subplots(figsize=figsize)

    cm = pd.DataFrame(cm)

    if use_unique:
        cm = cm.loc[s1.astype('category').cat.categories, s2.astype('category').cat.categories]

    if normalize:
        sns.heatmap(data=cm.dropna(), annot=True, ax=ax, fmt='0.2f', cmap=cmap, linewidth=linewidth, linecolor=linecolor, **kwargs)
    else:
        sns.heatmap(data=cm, annot=True, ax=ax, fmt='g', cmap=cmap, linewidth=linewidth, linecolor=linecolor, **kwargs)

    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.set_xlabel(s2.name, fontsize=figsize[1]*3)
    ax.set_ylabel(s1.name, fontsize=figsize[1]*3)

    for _, spine in ax.spines.items():
        spine.set_visible(True)

    return fig

def plot_cmatrix(
    cmatrix: pd.DataFrame,
    metric: str = 'euclidean',
    method: str = 'ward',
    n_clusters: int = 10,
    color_thresh_scale: float = 0.3,
    figsize: tuple = (8,8),
    p: int = 30,
    metas: Union[list, None] = None,
    vmax: Union[float, None] = None,
    vmin: Union[float, None] = None,
    cbar_label: str = 'ARD-NMF \nMembership',
    cmap: Union[str, None] = None,
    plot_cluster_lines: bool = False,
    cf: Union[None, pd.DataFrame] = None,
    cf_label: str = "Deconvolution",
    cf_cmap: str = "tab20c",
    show_xlab: bool = True,
    cdict_force: Union[dict,None] = None,
    dendro_flag: bool = False,
    dendro_meta: Union[pd.Series,None] = None,
    dendro_cmap: Union[dict,None] = None,
    dendro_groups_df: Union[pd.DataFrame,None] = None,
    dendro_pie_labels: Union[dict,None] = None,
    H_mut: pd.DataFrame = None,
    mut_colors: Union[dict,None] = None,
    show_upper_triangle = False,
    input_clusters = None,    
):
    """
    Plot consensus matrix.
    -----------------------
    Args:
        * cmatrix: consensus matrix. This may be generated by calling:
            df, assign_p = consensus_cluster_ardnmf(filepath)
        * metric: distance metric
        * method: method of clustering
        * n_clusters: number of clusters for agglomerative clustering
        * color_thresh_scale: asthetic scale for coloring of dendrogram
        * figsize: figsize
        * p: parameter for dendrogram
        * meta: list of pd.Series that includes a variable of interest to plot
            to left of plot; must be categorical in nature
    Returns:
        * fig
    """
    from scipy.cluster import hierarchy
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering
    import matplotlib as mpl
    from matplotlib.pyplot import cm
    import scipy.spatial.distance as ssd

    # -------------
    # Heatmap
    # -------------
    fig,ax = plt.subplots(figsize=figsize)
    cbar_ax = fig.add_axes([ax.get_position().x1 + ax.get_position().x1*0.01, ax.get_position().y0, .025, .1])

    # Compute initial linkage to grab ordering

    if metric=='precomputed':
        method='average'
        d_linkage =shc.linkage(shc.distance.squareform(1-cmatrix.values), method=method)
    else:
        d_linkage = shc.linkage(cmatrix, metric=metric, method=method)

    dres = shc.dendrogram(d_linkage, p=p, no_plot=True)
    dgram_idx = list(map(int, dres['ivl']))

    # Create heatmap
    if vmax is None:
        cbar_top_lim = np.max(cmatrix.values)
    else:
        cbar_top_lim = vmax

    if vmin is None:
        cbar_bottom_lim = 0
    else:
        cbar_bottom_lim = vmin

    if show_upper_triangle:
        mask = np.zeros_like(cmatrix, dtype=np.bool)
        mask[np.tril_indices_from(mask)] = True
        mask[np.diag_indices_from(mask)] = False
    else:
        mask=None

    sns.heatmap(
        cmatrix.iloc[dgram_idx,dgram_idx].values,
        ax=ax,
        square=True,
        cbar_ax=cbar_ax,
        cbar_kws = {'ticks':[cbar_bottom_lim, cbar_top_lim]},
        rasterized=True,
        vmax=vmax,
        vmin=vmin,
        cmap=cmap,
        mask=mask
    )

    cbar_ax.set_frame_on(True)
    cbar_ax.set_ylabel(cbar_label, fontsize=10, rotation=90)

    if show_xlab:
        ax.set_xticklabels(cmatrix.iloc[:,dgram_idx].columns, fontsize=8)
    else:
        ax.set_xticks([])
    ax.set_yticks([])

    x0 = ax.get_position().x0
    x1 = ax.get_position().x1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1

    buf = y1*0.015

    # Agglomerative Clustering
    cluster = AgglomerativeClustering(
        n_clusters=n_clusters,
        affinity=metric,
        linkage=method
    )

    clusters = cluster.fit_predict(cmatrix.iloc[dgram_idx,dgram_idx])
    cluster_color_list, _ = sa.pl.series_to_colors(pd.Series(clusters))

    # -------------------------
    # Mutational Signature Axis
    # -------------------------
    _end_height = y1+buf

    if isinstance(H_mut, pd.DataFrame):
        # Prep H_mut for plotting

        MUTSIG_CDICT={
            'HRD':'black', #(0.85, 0.11, 0.38)
            'POLE+MSI':'black', # (1, 0.76, 0.03)
            'Tobacco Smoking':'black', # (0.36, 0.86, 0.95)
            'APOBEC':'black', # (0, 0.30, 0.25)
            'TMB':'black' # (0.35, 0.0, 0.76)
        }

        H_mut = H_mut.loc[[x for x in H_mut.index if x in cmatrix.index.to_list()]]
        for x in cmatrix.index:
            if x not in H_mut.index.to_list():
                H_mut.loc[x] = 0
        patient_order = cmatrix.iloc[dgram_idx,dgram_idx].index
        H_mut = H_mut.loc[[x for x in patient_order if x in H_mut.index]]

        # Retrieve figure positions
        x0 = ax.get_position().x0
        x1 = ax.get_position().x1
        y0 = ax.get_position().y0
        y1 = ax.get_position().y1
        buf = y1*0.01

        # Add axes for mutational signature plots
        _height = 2.5*buf
        cohort_ax = fig.add_axes([x0,y1+0.5*buf, x1-x0, _height])
        hrd_cbar_ax = fig.add_axes([x0,y1+3.5*buf, x1-x0, _height])
        hrd_cbar_leg_ax = fig.add_axes([x1+buf,y1+3.5*buf,1.5*buf, _height])
        mmrd_cbar_ax = fig.add_axes([x0,y1+6.5*buf, x1-x0, _height])
        mmrd_cbar_leg_ax = fig.add_axes([x1+buf,y1+6.5*buf,1.5*buf, _height])
        smoking_cbar_ax = fig.add_axes([x0,y1+9.5*buf, x1-x0, _height])
        smoking_cbar_leg_ax = fig.add_axes([x1+buf,y1+9.5*buf,1.5*buf, _height])
        apobec_cbar_ax = fig.add_axes([x0,y1+12.5*buf, x1-x0, _height])
        apobec_cbar_leg_ax = fig.add_axes([x1+buf,y1+12.5*buf,1.5*buf, _height])
        tmb_cbar_ax = fig.add_axes([x0,y1+15.5*buf, x1-x0, _height])
        tmb_cbar_leg_ax = fig.add_axes([x1+buf,y1+15.5*buf,1.5*buf, _height])

        _end_height = y1+15.5*buf+_height+0.5*buf

        ## Cohort
        cluster_color_list, _ = sa.pl.series_to_colors(dendro_meta.iloc[dgram_idx], cdict=CPTAC_CMAP)
        mat, cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list)
        sns.heatmap(
            mat,
            cmap=cmap,
            ax=cohort_ax,
            xticklabels=False,
            yticklabels=False,
            cbar=False,
            rasterized=True
        )
        cohort_ax.set_ylabel("Cohort", rotation=0, y=0.25, ha='right', fontsize=16)
        [cohort_ax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

        ## HRD
        color_labels_hrd = H_mut['HRD'].astype(int).sort_values().unique()
        rgb_values_hrd = sns.light_palette(MUTSIG_CDICT['HRD'], len(color_labels_hrd))
        color_map_hrd = dict(zip(color_labels_hrd,rgb_values_hrd))

        # Patient to color map
        hrd_color_s = H_mut['HRD'].astype(int).map(color_map_hrd)

        # Create clor to attribution map
        hrd2color_dict = pd.concat([H_mut['HRD'].astype(int),hrd_color_s],1)
        hrd2color_dict.columns = ['HRD','color']
        hrd2color_dict = hrd2color_dict.drop_duplicates('color')
        hrd2color_dict = hrd2color_dict.set_index('color').sort_values(by='HRD').to_dict()['HRD']
        hrd_mat, hrd_cmap = sa.pl.color_list_to_matrix_and_cmap(hrd_color_s, hrd2color_dict)

        # Prep heatmap for log normalization
        hrd_mat = (hrd_mat[0]+1).reshape((1,len(hrd_mat[0])))
        sns.heatmap(
            hrd_mat,
            cmap=hrd_cmap,
            ax=hrd_cbar_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=hrd_cbar_leg_ax,
            norm=mpl.colors.LogNorm()
        )
        hrd_cbar_ax.set_ylabel("HRD", rotation=0, y=0.25, ha='right', fontsize=16)
        [hrd_cbar_ax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

        ## MMRD
        mmrd_s = H_mut['POLE+MSI'] + H_mut['CpG>T + Short Indels'] + H_mut['MMRD']
        color_labels_mmrd = mmrd_s.astype(int).sort_values().unique()
        rgb_values_mmrd = sns.light_palette(MUTSIG_CDICT['POLE+MSI'], len(color_labels_mmrd))
        color_map_mmrd = dict(zip(color_labels_mmrd,rgb_values_mmrd))

        # Patient to color map
        mmrd_color_s = mmrd_s.astype(int).map(color_map_mmrd)

        # Create color to attribution map
        mmrd2color_dict = pd.concat([mmrd_s.astype(int),mmrd_color_s],1)
        mmrd2color_dict.columns = ['MMRD','color']
        mmrd2color_dict = mmrd2color_dict.drop_duplicates('color')
        mmrd2color_dict = mmrd2color_dict.set_index('color').sort_values(by='MMRD').to_dict()['MMRD']
        mmrd_mat, mmrd_cmap = sa.pl.color_list_to_matrix_and_cmap(mmrd_color_s, mmrd2color_dict)

        # Prep heatmap for log normalization
        mmrd_mat = (mmrd_mat[0]+1).reshape((1,len(mmrd_mat[0])))
        sns.heatmap(
            mmrd_mat,
            cmap=mmrd_cmap,
            ax=mmrd_cbar_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=mmrd_cbar_leg_ax,
            norm=mpl.colors.LogNorm()
        )
        mmrd_cbar_ax.set_ylabel("MMRD", rotation=0, y=0.25, ha='right', fontsize=16)
        [mmrd_cbar_ax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

        ## Tobacco
        color_labels_smoking = H_mut['Tobacco Smoking'].astype(int).sort_values().unique()
        rgb_values_smoking = sns.light_palette(MUTSIG_CDICT['Tobacco Smoking'], len(color_labels_smoking))
        color_map_smoking = dict(zip(color_labels_smoking,rgb_values_smoking))
        smoking_color_s = H_mut['Tobacco Smoking'].astype(int).map(color_map_smoking)
        smoking2color_dict = pd.concat([H_mut['Tobacco Smoking'].astype(int),smoking_color_s],1)
        smoking2color_dict.columns = ['Tobacco Smoking','color']
        smoking2color_dict = smoking2color_dict.drop_duplicates('color')
        smoking2color_dict = smoking2color_dict.set_index('color').sort_values(by='Tobacco Smoking').to_dict()['Tobacco Smoking']
        smoking_mat, smoking_cmap = sa.pl.color_list_to_matrix_and_cmap(smoking_color_s, smoking2color_dict)
        smoking_mat = (smoking_mat[0]+1).reshape((1,len(smoking_mat[0])))

        sns.heatmap(
            smoking_mat,
            cmap=smoking_cmap,
            ax=smoking_cbar_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=smoking_cbar_leg_ax,
            norm=mpl.colors.LogNorm()
        )
        smoking_cbar_ax.set_ylabel("Tobacco Smoking", rotation=0, y=0.25, ha='right', fontsize=16)
        [smoking_cbar_ax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

        ## APOBEC
        color_labels_apobec = H_mut['APOBEC'].astype(int).sort_values().unique()
        rgb_values_apobec = sns.light_palette(MUTSIG_CDICT['APOBEC'], len(color_labels_apobec))
        color_map_apobec = dict(zip(color_labels_apobec,rgb_values_apobec))
        apobec_color_s = H_mut['APOBEC'].astype(int).map(color_map_apobec)
        apobec2color_dict = pd.concat([H_mut['APOBEC'].astype(int),apobec_color_s],1)
        apobec2color_dict.columns = ['APOBEC','color']
        apobec2color_dict = apobec2color_dict.drop_duplicates('color')
        apobec2color_dict = apobec2color_dict.set_index('color').sort_values(by='APOBEC').to_dict()['APOBEC']
        apobec_mat, apobec_cmap = sa.pl.color_list_to_matrix_and_cmap(apobec_color_s, apobec2color_dict)
        apobec_mat = (apobec_mat[0]+1).reshape((1,len(apobec_mat[0])))

        sns.heatmap(
            apobec_mat,
            cmap=apobec_cmap,
            ax=apobec_cbar_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=apobec_cbar_leg_ax,
            norm=mpl.colors.LogNorm()
        )

        apobec_cbar_ax.set_ylabel('APOBEC', rotation=0, y=0.25, ha='right', fontsize=16)
        [apobec_cbar_ax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

        ## TMB
        color_labels_tmb = H_mut['TMB'].astype(int).sort_values().unique()

        # Color
        rgb_values_tmb = sns.light_palette(MUTSIG_CDICT['TMB'], len(color_labels_tmb))

        color_map_tmb = dict(zip(color_labels_tmb,rgb_values_tmb))
        tmb_color_s = H_mut['TMB'].astype(int).map(color_map_tmb)
        tmb2color_dict = pd.concat([H_mut['TMB'].astype(int),tmb_color_s],1)
        tmb2color_dict.columns = ['TMB','color']
        tmb2color_dict = tmb2color_dict.drop_duplicates('color')
        tmb2color_dict = tmb2color_dict.set_index('color').sort_values(by='TMB').to_dict()['TMB']
        tmb_mat, tmb_cmap = sa.pl.color_list_to_matrix_and_cmap(tmb_color_s, tmb2color_dict)
        tmb_mat = (tmb_mat[0]+1).reshape((1,len(tmb_mat[0])))

        sns.heatmap(
            tmb_mat,
            cmap=tmb_cmap,
            ax=tmb_cbar_ax,
            yticklabels=False,
            xticklabels=False,
            cbar=True,
            cbar_ax=tmb_cbar_leg_ax,
            norm=mpl.colors.LogNorm()
        )

        tmb_cbar_ax.set_ylabel('TMB', rotation=0, y=0.25, ha='right', fontsize=16)
        [tmb_cbar_ax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

    # -------------
    # Dendrogram
    #
    # Or Deconvolution Results
    # -------------

    if cf is None:
        if dendro_flag:
            dax = fig.add_axes([x0, _end_height, x1-x0, 0.2]) ##
            plot_dendro_figure(
                cmatrix,
                dendro_meta,
                dendro_cmap,
                dendro_groups_df,
                color_thresh_scale=color_thresh_scale,
                pie_label=dendro_pie_labels,
                ax=dax
            )
        else:
            dax = fig.add_axes([x0, y1+buf, x1-x0, 0.15])
            cmap = cm.rainbow(np.linspace(0, 1, 10))
            hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])

            dres = shc.dendrogram(
                d_linkage,
                p=p,
                ax=dax,
                above_threshold_color="grey",
                color_threshold=color_thresh_scale*max(d_linkage[:,2])
            )

            dax.set_xticks([])
            dax.set_yticks([])
            [dax.spines[x].set_visible(False) for x in ['top','right','bottom','left']]
    else:
        dax = fig.add_axes([x0, y1+buf, x1-x0, 0.2])

        _d_idx = cmatrix.iloc[dgram_idx,dgram_idx].index
        _ = plot_cell_fraction(cf, order=_d_idx, vlines=np.unique(clusters,return_index=True)[1], ax=dax, cmap=cf_cmap)
        dax.set_xticks([])
        dax.set_yticks([])

        dax.set_ylabel(cf_label, fontsize=14)
        dax.set_xlabel("")


    # -------------
    # Metadata Axes
    # -------------
    if plot_cluster_lines:
        if input_clusters is None:
            input_clusters = clusters
        ###
        # plot horizontal lines
        hz_lines = np.sort(np.unique(pd.Series(input_clusters), return_index=True)[1])
        v,c = np.unique(input_clusters, return_counts=True)

        _c = hz_lines
        _c = np.roll(hz_lines, 1)
        _c[0] = 0
        _c[1] = 0

        _ci = hz_lines[1:]
        _ci = np.append(_ci, input_clusters.shape[0])

        if show_upper_triangle:
            for idx, hz in enumerate(hz_lines):
                ax.hlines(hz, hz, _ci[idx], rasterized=True)
                ax.vlines(_ci[idx], hz, _ci[idx], rasterized=True)
                ax.plot([1, 0], [0, 1], transform=ax.transAxes, c='black')
                #ax.plot([1, 0], [1, 1], transform=ax.transAxes, c='black')
                #ax.plot([1, 1], [1, 0], transform=ax.transAxes, c='black')
        else:
            for idx, hz in enumerate(hz_lines):
                ax.hlines(hz, _c[idx], _ci[idx], rasterized=True, color='yellow')
                ax.vlines(hz, _c[idx], _ci[idx], rasterized=True, color='yellow')


    # Add axes
    # Plots agglomerative clustering results
    if metas is None:
        lax = fig.add_axes([x0-3*buf, y0, 2*buf, y1-y0])
        mat, cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list)
        sns.heatmap(mat.T, cmap=cmap, ax=lax, xticklabels=False, yticklabels=False, cbar=False, rasterized=True)

        uniq, idx, num_vals = np.unique(clusters.T, return_index=True, return_counts=True)
        y_locs = idx + num_vals / 2

        for idx,u in enumerate(uniq):
            lax.text(x0-50*buf, y_locs[idx], u, ha='center')

        for idx,u in enumerate(uniq):
            ax.text(
                mat.shape[1]+0.01*mat.shape[1],
                y_locs[idx],
                "n={}".format(num_vals[idx]),
                ha='left',
                fontsize=14
            )

        for _, spine in lax.spines.items():
            spine.set_visible(True)

        lax.set_xlabel("Consensus", rotation=90)

    else:
        for idx,meta in enumerate(metas):
            new_ax = [x0-(idx+3)*buf-(idx*2)*buf, y0, 2*buf, y1-y0]
            lax = fig.add_axes(new_ax)

            if isinstance(meta, str) and meta=='aggr':
                mat, cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list)
                sns.heatmap(mat.T, cmap=cmap, ax=lax, xticklabels=False, yticklabels=False, cbar=False, rasterized=True)

                uniq, idx, num_vals = np.unique(clusters.T, return_index=True, return_counts=True)
                y_locs = idx + num_vals / 2

                for idx,u in enumerate(uniq):
                    lax.text(0.5, y_locs[idx], u, ha='center')

                for idx,u in enumerate(uniq):
                    ax.text(
                        mat.shape[1]+0.01*mat.shape[1],
                        y_locs[idx],
                        "n={}".format(num_vals[idx]),
                        ha='left',
                        fontsize=14
                    )

                lax.set_xlabel("cluster", rotation=90)

            else:
                if meta.name == 'cohort':
                    cdict = CPTAC_CMAP
                    order_dict = None
                    #order_dict = {COLOR_SCHEME['pink']:0, COLOR_SCHEME['yellow']:1, COLOR_SCHEME['teal']:2, COLOR_SCHEME['purple']:3,  COLOR_SCHEME['orange']:4,  COLOR_SCHEME['lightblue']:5}
                elif meta.name == 'staging':
                    cdict = TUMOR_STAGE_SCHEME
                    order_dict = None
                else:
                    if cdict_force is not None:
                        cdict = cdict_force
                    else:
                        cdict = None
                    order_dict = None

                meta = meta.loc[cmatrix.index[dgram_idx]]
                _meta_v = pd.unique(meta)

                if None in _meta_v and True in _meta_v and len(_meta_v)==2:
                    cdict = {True:(0.21044753832183283, 0.6773105080456748,0.6433941168468681), None:(1.0,1.0,1.0,1.0)}
                    order_dict = {(1.0,1.0,1.0,1.0):0,(0.21044753832183283,0.6773105080456748,0.6433941168468681):1}

                cluster_color_list, _ = sa.pl.series_to_colors(meta,cdict=cdict)
                mat,cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list,order_dict=order_dict)
                sns.heatmap(mat.T, cmap=cmap, ax=lax, yticklabels=False, xticklabels=False, cbar=False)
                lax.set_xlabel(meta.name, rotation=90)

            for _, spine in lax.spines.items():
                spine.set_visible(True)

    rs = pd.DataFrame(clusters, index=cmatrix.index[dgram_idx]).rename(columns={0:'clusters'})

    for _side, spine in ax.spines.items():
        if show_upper_triangle:
            if _side in ['Top','Right']:
                spine.set_visible(True)
        else:
            spine.set_visible(True)

    if show_xlab:
        ax.set_xlabel("Samples", fontsize=14)
    else:
        ax.set_xlabel("")



    return fig, rs, dres, d_linkage

def drawPieMarker(xs, ys, ratios, sizes, colors, ax=None, zorder=1, resolution=40):
    """
    Draw Pie Marker

    """
    if ax is None:
        ax = plt.gca()

    markers = []
    previous = 0
    for color, ratio in zip(colors, ratios):
        this = 2 * np.pi * ratio + previous
        x  = [0] + np.cos(np.linspace(previous, this, resolution)).tolist() + [0]
        y  = [0] + np.sin(np.linspace(previous, this, resolution)).tolist() + [0]
        xy = np.column_stack([x, y])
        previous = this
        markers.append({'marker':xy, 's':np.abs(xy).max()**2*np.array(sizes), 'facecolor':color, 'edgecolor':'k'})

    # scatter each of the pie pieces to create pies
    for marker in markers:
        ax.scatter(xs, ys, zorder=zorder, **marker)

def plot_dendro_figure(
    cmatrix: pd.DataFrame,
    pie_meta: pd.Series,
    pie_cmap: dict,
    groups_df: pd.DataFrame,
    metric: str = 'euclidean',
    method: str = 'ward',
    figsize: tuple = (20,6),
    color_thresh_scale: float = 0.3,
    p: int = 30,
    pie_chart_size: int = 250,
    pie_label: Union[None, str, dict] = "id",
    adjust_txt: bool = False,
    ax: Union[None, plt.axes] = None,
    pie_label_small: Union[None, list] = None
    ):
    """
    Plot dendrogram figure.
    """
    from scipy.cluster import hierarchy
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering
    import matplotlib as mpl
    from matplotlib.pyplot import cm
    import scipy.spatial.distance as ssd
    from adjustText import adjust_text

    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    d_linkage = shc.linkage(cmatrix, metric=metric, method=method)
    dres = shc.dendrogram(
        d_linkage,
        p=p,
        ax=ax,
        above_threshold_color="grey",
        color_threshold=color_thresh_scale*max(d_linkage[:,2])
    )

    ax.set_xticks([])
    ax.set_yticks([])
    [ax.spines[x].set_visible(False) for x in ['top','right','bottom','left']]

    dgram_idx = list(map(int, dres['ivl']))
    A = np.array(dres['dcoord'])
    B = np.array(dres['icoord'])
    _order = A[:, 1].argsort()[::-1]

    c = 1
    texts = list()
    for i, d in zip(B[_order], A[_order]):
        if c in groups_df:
            x = 0.5 * sum(i[1:3])
            y = d[1]
            l,n = np.unique(groups_df.loc[:,[c]].dropna().join(pie_meta)[pie_meta.name],return_counts=True)
            drawPieMarker(x, y, n / np.sum(n), sizes=[pie_chart_size], colors=[pie_cmap[x] for x in l], ax=ax, zorder=2)

            if pie_label is not None:
                if pie_label=="id":
                    _text = "{}".format(c)
                else:
                    _text = "{}".format(pie_label[c])

                if pie_label_small is not None:
                    if c in [9,12]:
                        t = ax.text(x-.5, y+6, _text, ha='center', fontsize=8)
                    else:
                        t = ax.text(x-.5, y+6, _text, ha='center', fontsize=14)
                else:
                    t = ax.text(x-.5, y+6, _text, ha='center', fontsize=8)

                t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

                texts.append(t)

        c+=1

    if adjust_txt:
        #pass
        adjust_text(texts, ax=ax, only_move={'points':'x', 'text':'x'})

    ax.set_xlim(*ax.get_xlim())
    ax.set_ylim(*ax.get_ylim())

    for l,k in pie_cmap.items():
        ax.scatter(-5,-5, c=np.array([k]), label=l)

    ax.legend(frameon=False, ncol=2, loc='center left', bbox_to_anchor=(1, 0.5))

def plot_enrichment(
    df,
    ax=None,
    padj_thresh=1e-2,
    figsize=(6,6),
    size_min=5,
    labelsize=9,
    cmap='cool',
    title=None,
    padj_idx='padj',
    cbar_ax_title='Adj P-value',
    sort_by=None,
    sort_by_ascending=True,
    xlim=None,
    xax='ES',
    n_to_plot=15,
    size_multi=1.25
    ):
    """
    Plot GSEA results
    -----------------------
    Args:
        * df: pd.DataFrame
        * ax: matplotli.pyplot.axis
        * padj_thresh: threshold value to plot for `padj`
        * figsize: tuple
        * size_min: minimum gene overlap in pathway
        * labelsize: int size of pathway labels
        * cmap: colormap

    Returns:
        * fig
    """
    from matplotlib import rcParams
    from matplotlib import colors

    if ax == None:
        fig,ax = plt.subplots(figsize=figsize)

    if sort_by is None:
        sort_by = padj_idx
        sort_by_ascending = True

    data_to_plot = df[(df[padj_idx] < padj_thresh) & (df['size']>size_min)].copy()
    data_to_plot = data_to_plot.sort_values(sort_by, ascending=sort_by_ascending)

    if sort_by == 'padj':
        # Plot top by significance
        data_to_plot = data_to_plot.head(n_to_plot*2)
    elif sort_by in ("NES","ES"):
        # Plot top / bottom by ES
        data_to_plot = pd.concat((data_to_plot.head(n_to_plot), data_to_plot.tail(n_to_plot))).drop_duplicates(subset='pathway')

    xlim_bound = max(np.max(data_to_plot[xax]), np.abs(np.min(data_to_plot[xax])))

    norm = colors.LogNorm(data_to_plot[padj_idx].min(), data_to_plot[padj_idx].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    rcParams.update({'font.size': 14, 'font.weight': 'normal'})

    ax.set_axisbelow(True)

    if data_to_plot.shape[0]==0:
        print("   * no significant pathways")
        return
    else:
        print("   * plotting {} pathways".format(data_to_plot.shape[0]))

    path = ax.scatter(
        x=xax,
        y='pathway',
        c=padj_idx,
        cmap=cmap,
        norm=norm,
        data=data_to_plot,
        linewidth=1,
        edgecolor="black",
        s=[(i+10)**size_multi for i in data_to_plot['size']],
        zorder=10
    )

    # Plot Bars
    bar_colors = ['pink' if x > 0 else 'lightblue' for x in data_to_plot[xax]]
    data_to_plot.plot(
        kind='barh',
        y=xax,
        x='pathway',
        legend=None,
        ax=ax,
        alpha=0.5,
        color=bar_colors,
        width=0.7,
        zorder=1,
    )


    # Lower p-value on top
    ax.invert_yaxis()
    ax.tick_params(axis='y', which='major', labelsize=labelsize)
    ax.set_ylabel('')
    ax.set_xlabel(xax, fontsize=14, fontweight='normal')
    ax.xaxis.grid(False)
    ax.yaxis.grid(True)

    if xlim is None:
        ax.set_xlim(-xlim_bound-.1, xlim_bound+.1)
    else:
        ax.set_xlim(xlim)

    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # -------------------------
    # p-value Colorbar
    # -------------------------
    cbaxes = ax.figure.add_axes([0.8, 0.125, 0.03, 0.2])
    cbar = ax.figure.colorbar(sm, shrink=0.4, cax=cbaxes)
    cbar.ax.invert_yaxis()

    cbar.ax.set_title(cbar_ax_title, fontweight='normal', ha='center', x=2, y=1.15, fontsize=12)
    cbar.ax.tick_params(axis='y', which='major', labelsize=8)
    cbar.ax.tick_params(axis='y', which='minor', labelsize=8)

    # -------------------------
    # Size Legend
    # -------------------------
    min_bound = int(np.floor(min(data_to_plot['size']) / 10.0)) * 10
    max_bound = int(np.ceil(max(data_to_plot['size']) / 10.0)) * 10
    med_bound = int((min_bound + max_bound) / 2)

    l1 = ax.scatter([],[], s=(min_bound+10)**size_multi, edgecolors='none', color='black')
    l2 = ax.scatter([],[], s=(med_bound+10)**size_multi, edgecolors='none', color='black')
    l3 = ax.scatter([],[], s=(max_bound+10)**size_multi, edgecolors='none', color='black')
    labels = [min_bound, med_bound, max_bound]

    leg = ax.legend(
        [l1, l2, l3],
        labels,
        ncol=1,
        frameon=False,
        fontsize=12,
        handlelength=2.5,
        loc = 'upper right',
        labelspacing = 2,
        bbox_to_anchor=(1.45,1),
        title='Gene Overlap',
        scatterpoints = 1,
        facecolor='black'
    )

    ax.axvline(0, color='black', linewidth=1, alpha=0.5)
    ax.set_title(title)

    return data_to_plot

def plot_ptm_de_proportions(
    up_df,
    down_df,
    figsize=(16,4)
    ):
    """
    Plot PTM DE Proportions
    -------------------------
    Args:
        * up_df: pd.DataFrame
        * down_df: pd.DataFrame
        * figsize: tuple of figsize

    Returns:
        * fig
    """
    fig,axes = plt.subplots(1, 2, figsize=figsize)

    plot_proportions(up_df, color=[PTM_SCHEME2[x.capitalize()] for x in up_df.columns], figsize=(4,4), ax=axes[0])
    plot_proportions(down_df, color=[PTM_SCHEME2[x.capitalize()] for x in down_df.columns], figsize=(4,4), ax=axes[1])

    axes[0].legend().remove()
    axes[0].set_xlabel('Up', fontsize=12)
    axes[0].set_ylabel(up_df.index.name, fontsize=12)
    axes[1].set_xlabel('Down', fontsize=12)
    axes[1].set_ylabel('')

    return fig

def _add_dendrogram(
    X,
    dax=None,
    metric: str = 'euclidean',
    method: str = 'ward',
    p: int = 30,
    color_thresh_scale: float = 0,
    ):
    """
    Compute feature order and plot dendrogram.
    """
    from scipy.cluster import hierarchy as shc
    from scipy.cluster import hierarchy
    from sklearn.cluster import AgglomerativeClustering
    import matplotlib as mpl
    from matplotlib.pyplot import cm

    d_linkage = shc.linkage(X, metric=metric, method=method)
    dres = shc.dendrogram(d_linkage, p=p, no_plot=True)
    dgram_idx = list(map(int, dres['ivl']))

    if dax is None:
        return dgram_idx

    cmap = cm.rainbow(np.linspace(0, 1, 10))
    hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])

    dres = shc.dendrogram(
        d_linkage,
        p=p,
        ax=dax,
        above_threshold_color="grey",
        color_threshold=color_thresh_scale*max(d_linkage[:,2]),
        orientation='left'
    )

    dax.set_xticks([])
    dax.set_yticks([])
    [dax.spines[x].set_visible(False) for x in ['top','right','bottom','left']]

    return dgram_idx

def _add_metadata_bar(meta, ax, cdict=None, order_dict=None, label=None):
    """
    Add metadata barplot. CPTAC.
    """
    from plotting import COLOR_SCHEME, TUMOR_STAGE_SCHEME, CPTAC_CMAP, PTM_SCHEME2
    import signatureanalyzer as sa

    if meta.name == 'Feature':
        cdict = PTM_SCHEME2

    cluster_color_list, _ = sa.pl.series_to_colors(meta, cdict=cdict)
    mat, cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list, order_dict=order_dict)
    sns.heatmap(mat, cmap=cmap, ax=ax, yticklabels=False, xticklabels=False, cbar=False)

    if label is None:
        label = meta.name

    ax.yaxis.set_label_position("right")
    ax.set_ylabel(label, rotation=360, ha='left', va='center')
    ax.set_yticks([])

def plot_cptac_pathways_heatmap(
    e_df,
    cluster_idx='id',
    feature_idx='feature',
    pathway_idx='pathway',
    p_adj_idx='padj',
    es_idx='NES',
    clust_order=None,
    feat_order=None,
    by_cluster=False,
    figsize=(12,36),
    padj_thresh=0.1,
    show_pval=True
    ):
    """
    Plot CPTAC pathways heatmap.
    """
    from matplotlib.gridspec import GridSpec
    from matplotlib.patches import Rectangle

    if clust_order is None:
        clust_order = [3,5,4,2,1,0]
    if feat_order is None:
        feat_order = ['Acetylome','Phosphoproteome','Proteome','Transcriptome']

    e_df['lab'] = e_df[cluster_idx].astype(str) + "-" + e_df[feature_idx]

    # Sample indices
    if by_cluster:
        samp_idx = [str(c)+'-'+f for c in clust_order for f in feat_order]
    else:
        samp_idx = [str(c)+'-'+f for f in feat_order for c in clust_order]

    si_df = pd.DataFrame(samp_idx, columns=['id'])
    si_df['Cluster'] = si_df['id'].apply(lambda x: x.split('-')[0])
    si_df['Feature'] = si_df['id'].apply(lambda x: x.split('-')[1])
    si_df = si_df.set_index('id')

    # Create dataframes
    padj_df = e_df.pivot(index=pathway_idx, values=p_adj_idx, columns='lab').loc[:,si_df.index]
    nes_df = e_df.pivot(index=pathway_idx, values=es_idx, columns='lab').loc[:,si_df.index]

    # Filter
    mask = padj_df.min(1)<padj_thresh
    nes_df = nes_df.loc[mask]
    padj_df = padj_df.loc[mask]

    fig = plt.figure(figsize=figsize)

    gs = GridSpec(nes_df.shape[0]+1,12, figure=fig)
    gs.update(wspace=0.15, hspace=0.15)

    dax = fig.add_subplot(gs[2:,:2])
    ax = fig.add_subplot(gs[2:,2:])

    dgram_idx = _add_dendrogram(nes_df.fillna(0), dax)
    sample_idx = _add_dendrogram(nes_df.fillna(0).T)

    # Add colorbar
    cbar_ax = fig.add_axes([1.25, 0.13, .025, .1])

    # --------------------------------
    # Heatmap
    # --------------------------------
    sns.heatmap(
        nes_df.iloc[dgram_idx],
        ax=ax,
        cmap='coolwarm',
        cbar_ax=cbar_ax
    )

    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.yaxis.tick_right()
    _ = ax.set_yticks(np.arange(0,len(dgram_idx))+0.5)
    _ = ax.set_yticklabels(nes_df.index[dgram_idx], rotation=360)

    ax.xaxis.tick_top()
    n_ticks = np.sort(np.unique(si_df['Cluster'], return_index=True)[1])+2
    n_labs = pd.unique(si_df['Cluster'])

    ax.set_xticks([])


    feat_ax = fig.add_subplot(gs[0,2:])
    clust_ax = fig.add_subplot(gs[1,2:])

    _add_metadata_bar(si_df['Cluster'], clust_ax)
    _add_metadata_bar(si_df['Feature'], feat_ax)

    if show_pval:
        x,y = np.where(padj_df.iloc[dgram_idx].values<=0.1)

        print(x.shape,y.shape)

        for xi,yi in zip(x,y):
            _adj_p = -np.log10(padj_df.iloc[dgram_idx].iloc[xi,yi])
            ax.add_patch(Rectangle((yi,xi), 1, 1, fill=False, edgecolor='yellow', lw=_adj_p/2))

    cbar_ax.set_title(es_idx)
    return nes_df.iloc[dgram_idx],padj_df.iloc[dgram_idx]

def plot_volcano(
    de,
    thresh=20,
    ax=None,
    xlim=None,
    yax='qval',
    xax='logFC',
    xlabel=None,
    ylabel=None,
    gene_id='genes',
    ptm_id=None,
    filter_noise=False,
    fix_extremes=True,
    arrow=False,
    label_percentile=99.5,
    shuffle=True,
    gene_fontsize=10,
    only_plot_gene_symbol=True,
    genes_to_label=None,
    c1='black',
    c2='red'
    ):
    """
    Differential Expression Volcano Plot
    ---------------------
    Args:
        * de: pd.DataFrame of differentail expression results
        * thresh: qval threshold for differential expression
        * ax: matplotlib.axis
        * xlim: x-limits for the volcano plot; automatically determined / centered
            if not-provided
        * yax: column in dataframe for y-axis ('qval')
        * xax: column in dataframe for x-axis (usually log-fold-change)
        * xlabel: plot x-label
        * ylabel: plot y-label
        * gene_id: column in dataframe to use for gene-labels
        * ptm_id: column in dataframe to use for assay type if necessary
        * filter_noise: remove noise (high LFC, zero significance)
        * fix_extremes: fix-extremes
        * arrow: add arrows to each gene label leading back to point
        * label_percentile: percentile of genes on each neg/pos to label
        * shuffle: shuffle points (relevant for different assays on same plot)
        * gene_fontsize: fontsize for each gene
        * only_plot_gene_symbol: plot gene symbol (if assays available, will plot PTM modification)

    Returns:
        * returns fig
    """
    from numpy import inf
    from adjustText import adjust_text

    # Subset relevant Columns
    if ptm_id in de:
        de = de.loc[:,[xax, yax, gene_id, ptm_id, 'variableSites']]
    else:
        de = de.loc[:,[xax, yax, gene_id]]

    if xlim is None:
        _x = np.ceil(np.max(np.abs(de[xax])))
        xlim = (-_x, _x)

    # Shuffle columns
    if shuffle:
        de = de.sample(frac=1)

    PTM_VP_SCHEME = {
        'Proteome': (0.9677975592919913, 0.44127456009157356, 0.5358103155058701),
        'Phosphoproteome': (0.3126890019504329, 0.6928754610296064, 0.1923704830330379),
        'Acetylome': (0.23299120924703914, 0.639586552066035, 0.9260706093977744)
    }

    # Create -log10q-value
    de['-logQ'] = -np.log10(de[yax])

    # Get percentile threshold for each
    de_up = de[de[xax]>0]
    de_down = de[de[xax]<0]

    # Filter extremes
    if fix_extremes:
        def _fix(x,v):
            if x==inf:
                return int(v)
            else:
                return x

        val = int(np.max(de[de['-logQ'] != inf]['-logQ'])+1)
        de['-logQ'] = de['-logQ'].apply(lambda x: _fix(x,val))

    if filter_noise:
        de[(de[xax] < -20) & (de[yax] < 1e-3)]

        noisy_down_genes = de[(de[xax]<-5) & (de[yax]>1e-3)][gene_id]
        noisy_up_genes = de[(de[xax]>5) & (de[yax]>1e-3)][gene_id]
        noisy_genes = set(noisy_up_genes) | set(noisy_down_genes)

        de = de[~de[gene_id].isin(noisy_genes)]

    lowqval_de = de.loc[de['-logQ'] > thresh]
    other_de = de.loc[de['-logQ'] < thresh]

    if ax is None:
        fig, ax = plt.subplots(figsize=(6,6))

    # Below Threshold
    sns.regplot(
        other_de[xax],
        np.abs(other_de['-logQ']),
        fit_reg=False,
        scatter_kws={'s':12, 'alpha':0.5, 'color':c1,'rasterized':True},
        ax=ax,
    )

    # Above Threshold
    if ptm_id is None:
        sns.regplot(
            lowqval_de[xax],
            np.abs(lowqval_de['-logQ']),
            fit_reg=False,
            scatter_kws={'s':12, 'alpha':0.5,'color':c2,'rasterized':True},
            ax=ax,
        )
    else:
        # Split by PTM
        for ptm in pd.unique(lowqval_de[ptm_id]):
            _lowqval_de_ptm = lowqval_de[lowqval_de[ptm_id]==ptm]
            sns.regplot(
                _lowqval_de_ptm[xax],
                np.abs(_lowqval_de_ptm['-logQ']),
                fit_reg=False,
                scatter_kws={'s':12, 'alpha':0.5,'color':PTM_VP_SCHEME[ptm]},
                ax=ax,
            )

    # Axis labels
    if xax=='diff':
        xax_label = "Difference in Means"
    else:
        xax_label = xax.capitalize()

    ax.set_xlabel(xax_label, fontsize=14)
    ax.set_ylabel("-log10 {}".format(yax.capitalize()), fontsize=14)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))

    # Label names and positions
    #lowqval_de = lowqval_de.dropna()

    if ptm_id is None:
        # No PTM
        _up_thresh = np.percentile(de_up['-logQ'], label_percentile)
        _down_thresh = np.percentile(de_down['-logQ'], label_percentile)

        x = lowqval_de[xax].values
        y = lowqval_de['-logQ'].values
        labels = lowqval_de[gene_id].values

        to_keep = []
        for idx,x_val in enumerate(x):
            if x[idx] < 0 and y[idx] > _down_thresh:
                to_keep.append(idx)
            elif x[idx] > 0 and y[idx] > _up_thresh:
                to_keep.append(idx)

        x = x[to_keep]
        y = y[to_keep]
        labels = labels[to_keep]
    else:
        # Within each PTM
        x = list()
        y = list()
        labels = list()

        for _ptm in pd.unique(lowqval_de[ptm_id]):
            _lowqval_de_ptm = lowqval_de[lowqval_de[ptm_id]==_ptm]

            _up_thresh = np.percentile(de_up[de_up[ptm_id]==_ptm]['-logQ'], label_percentile)
            _down_thresh = np.percentile(de_down[de_down[ptm_id]==_ptm]['-logQ'], label_percentile)

            _x = _lowqval_de_ptm[xax].values
            _y = _lowqval_de_ptm['-logQ'].values

            if only_plot_gene_symbol:
                _labels = _lowqval_de_ptm[gene_id].values
            else:
                if _ptm=='Proteome':
                    _labels = _lowqval_de_ptm[gene_id].values
                else:
                    _labels = np.array(_lowqval_de_ptm[gene_id] + "\n("+_lowqval_de_ptm.variableSites.str.strip().fillna("")+")")

            to_keep = []
            for idx,x_val in enumerate(_x):
                if _x[idx] < 0 and _y[idx] > _down_thresh:
                    to_keep.append(idx)
                elif _x[idx] > 0 and _y[idx] > _up_thresh:
                    to_keep.append(idx)

            x += list(_x[to_keep])
            y += list(_y[to_keep])
            labels += list(_labels[to_keep])

    # Additional genes to label
    if genes_to_label is not None:
        _de_to_label = de[de[gene_id].isin(genes_to_label)]
        x = np.concatenate((x,_de_to_label[xax].values))
        y = np.concatenate((y,_de_to_label['-logQ'].values))
        labels = np.concatenate((labels,_de_to_label[gene_id].values))

    ax.axvline(0, color='black', alpha=0.2, rasterized=True)
    ax.axhline(0, color='black', alpha=0.2, rasterized=True)

    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=14)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=14)

    # Set xlabel
    if xlim is not None:
        ax.set_xlim(xlim)
    else:
        xlim = ax.get_xlim()

    texts = list()
    for i,txt in enumerate(labels):
        if x[i] > xlim[0] and x[i] < xlim[1]:
            texts.append(ax.text(x[i], y[i], txt, ha='center', va='center', fontsize=gene_fontsize))
    if arrow:
        adjust_text(texts, arrowprops=dict(arrowstyle='->', color='blue'))
    else:
        adjust_text(texts)

    return ax

def plot_de_qq_volcano(
    df,
    groupby='consensus',
    save_dir=None,
    pval_idx='pval',
    **kwargs
    ):
    """
    Plot.

    Plot DE QQ Volcano Args:
        * df: differential expression dataframe
        * groupby: groupings
        * save_dir: directory to save results in
        * pval_idx: for qq plot

    kkwargs for plot_volcano:
        * de: pd.DataFrame of differentail expression results
        * thresh: qval threshold for differential expression
        * ax: matplotlib.axis
        * xlim: x-limits for the volcano plot; automatically determined / centered
            if not-provided
        * yax: column in dataframe for y-axis ('qval')
        * xax: column in dataframe for x-axis (usually log-fold-change)
        * xlabel: plot x-label
        * ylabel: plot y-label
        * gene_id: column in dataframe to use for gene-labels
        * ptm_id: column in dataframe to use for assay type if necessary
        * filter_noise: remove noise (high LFC, zero significance)
        * fix_extremes: fix-extremes
        * arrow: add arrows to each gene label leading back to point
        * label_percentile: percentile of genes on each neg/pos to label
        * shuffle: shuffle points (relevant for different assays on same plot)
        * gene_fontsize: fontsize for each gene
        * only_plot_gene_symbol: plot gene symbol (if assays available, will plot PTM modification)
    """
    groups = pd.unique(df[groupby])

    for idx,group in enumerate(groups):
        fig,axes = plt.subplots(1, 2, figsize=(12, 6))
        df_group = df[df[groupby]==group]

        _ = plot_volcano(
            df_group,
            ax=axes[0],
            **kwargs
        )

        qqplot(df_group[pval_idx], ax=axes[1])
        plt.suptitle("Cluster {}".format(group), fontsize=16, x=0.125, y=1.025, fontweight="bold")
        plt.tight_layout()

        if save_dir is not None:
            os.makedirs(save_dir, exist_ok=True)
            plt.savefig(os.path.join(save_dir, "{}_{}_de.pdf".format(groupby, group)), dpi=200, bbox_inches='tight')

def plot_sankey(df, var1='phi1', var2='phi2', count='n', outfile=None, return_gb=False, **kwargs):
    """
    Generate Sankey Diagram.
    --------------------------
    Args:
        * df: pd.DataFrame
        * var1: column of df
        * var2: colmnn of df
            Creates relationship from var1->var2
        * count: name given to count column of dataframe
        * outfile: filename to save figure
        * return_gb: returns grouped dataframe with index(samples) for each gruoping
        * kkwargs: for creating node with go.Sankey
    """
    import plotly.graph_objects as go

    assert var1 in df, "{} not in input dataframe.".format(var1)
    assert var2 in df, "{} not in input dataframe.".format(var2)

    gb_df = df.groupby([var1,var2]).size().reset_index().rename(columns={0:count})

    sankey_map1 = {y:x for x,y in enumerate(np.unique(gb_df[var1]))}
    sankey_map2 = {y:x+np.max(list(sankey_map1.values()))+1 for x,y in enumerate(np.unique(gb_df[var2]))}

    gb_df[str(var1)+'_sankey'] = gb_df.loc[:,var1].apply(lambda x: sankey_map1[x])
    gb_df[str(var2)+'_sankey'] = gb_df.loc[:,var2].apply(lambda x: sankey_map2[x])

    labels = [str(x) for x in np.unique(gb_df[var1])] + [str(x) for x in np.unique(gb_df[var2])]

    link = dict(
        source=gb_df.loc[:,str(var1)+'_sankey'],
        target=gb_df.loc[:,str(var2)+'_sankey'],
        value=gb_df.loc[:,count]
    )

    node = dict(label=labels, **kwargs)

    data = go.Sankey(link = link, node=node)

    # plot
    fig = go.Figure(data)
    fig.show()

    if outfile is not None:
        fig.write_image(outfile)

    if return_gb:
        gb_full_df = gb_df[[var1,var2,count]]

        for idx,row in gb_full_df.iterrows():
            _l = list(df[(df[var1]==row[var1]) & (df[var2]==row[var2])].index)
            gb_full_df.loc[idx,'samples'] = ','.join(_l)

        return gb_full_df

def plot_signatures_heatmap(nmf_file, metas, order=None, cmap='Blues', figsize=(8,8)):
    """
    Plot signatures heatmap.
    """
    # Get ordering
    if order is None:
        order = pd.read_hdf(nmf_file, "H").sort_values(["max_id","max_norm"], ascending=[True,False]).index

    # Load H-matrix
    H = pd.read_hdf(nmf_file, "Hraw").loc[:,order]

    # Get colors
    meta_series_colors = list()
    for meta in metas:
        if meta.name == 'cohort':
            meta_series_colors.append(meta.loc[H.columns].apply(lambda x: CPTAC_CMAP[x]))
        elif meta.name == 'staging':
            meta_series_colors.append(meta.loc[H.columns].apply(lambda x: TUMOR_STAGE_SCHEME[x]))

    # Clustermap
    g = sns.clustermap(
        H,
        col_colors=meta_series_colors,
        figsize=figsize,
        cmap=cmap,
        cbar_kws={"label":"Signature Weight",},
        cbar_pos=(1.025, 0.175, .03, .2),
        col_cluster=False,
        row_cluster=False,
        yticklabels=True
    )

    ax = g.ax_heatmap
    ax.set_xticks([])
    ax.set_xlabel('')

    [ax.spines[s].set_visible(True) for s in ['top','bottom','left','right']]

#     for label in metas[0].unique():
#         if not isinstance(label, float):
#             g.ax_col_dendrogram.bar(0, 0, color=pl.TUMOR_STAGE_SCHEME[label], label=label, linewidth=0)
#     g.ax_col_dendrogram.legend(bbox_to_anchor=(-0.025, 0.05), loc="left", ncol=1)

    for label in metas[0].unique():
        if not isinstance(label, float):
            g.ax_col_dendrogram.bar(0, 0, color=CPTAC_CMAP[label], label=label, linewidth=0)
    g.ax_col_dendrogram.legend(bbox_to_anchor=(-0.025, 0.05), loc="left", ncol=1, frameon=False)

def plot_strip(
    X: pd.DataFrame,
    group_s: pd.Series,
    yax: str,
    figsize: tuple = (10,3),
    s: float = 5,
    order: Union[None, list] = None,
    ylabel: Union[None, str] = None,
    plot_violin: bool = True,
    plot_mean_line: bool = True,
    y_offset: float = 0.5,
    ylim: Union[None, tuple] = None,
    show_pval: Union[str,None] = None,
    ax = None,
    format_fig: bool = False,
    show_legend: bool = False,
    annot = None,
    color = None,
    **kwargs
    ):
    """
    Plot Strip Plot.
    -----------------------
    Args:
        * X: pd.DataFrame (features x samples)
        * group_s: pd.Series with annotations for samples
        * yax: feature to plot
        * figsize: tuple for figure size
        * s: size of dots
        * order: ordering of subtypes
        * ylabel: str
    Returns:
        * figure
    """
    from scipy.stats import kruskal

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    categ = list(set(group_s))
    pos = [x+1 for x in range(len(categ))]

    split_by = group_s.name
    _X = X.T.join(group_s)

    kruskal_statistic, kruskal_pval = kruskal(*_X.groupby(split_by)[yax].apply(np.array).values)

    if color is None:
        color = sns.color_palette("husl", len(categ))

    sns.stripplot(
        x=split_by,
        y=yax,
        data=_X,
        ax=ax,
        s=s,
        palette=color,
        order=order,
        **kwargs
    )

    if plot_mean_line:
        # Add means
        capprops = dict(linewidth=0, color='white')
        sns.boxplot(
            x=split_by,
            y=yax,
            data=_X,
            ax=ax,
            showfliers=False,
            showbox=False,
            whis=0,
            capprops=capprops,
            order=order
        )

    if plot_violin:
        sns.violinplot(
            x=split_by,
            y=yax,
            data=_X,
            ax=ax,
            alpha=0.2,
            color='white',
            order=order
        )

    if ylim is None:
        ax.set_ylim((ax.get_ylim()[0]-y_offset,ax.get_ylim()[1]+y_offset))
    else:
        ax.set_ylim(ylim)

    # Y label
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    # Add axis
    leg = ax.legend([yax], handlelength=0, handletextpad=0, fancybox=True, fontsize=12)
    for item in leg.legendHandles:
        item.set_visible(False)

    # Add p-value
    if show_pval is not None:
        if show_pval == 'lower left':
            ax.annotate("   Kruskal-Wallis, p = {:.3E}\n".format(kruskal_pval),(ax.get_xlim()[0],ax.get_ylim()[0]),fontsize=8)
        elif show_pval == 'above left':
            ax.annotate("   Kruskal-Wallis, p = {:.3E}\n".format(kruskal_pval),(ax.get_xlim()[0],ax.get_ylim()[1]),fontsize=8)

    if format_fig:
        _format_ax(ax)

    if annot is not None:
        ax.annotate("{}\n".format(annot),(ax.get_xlim()[1], ax.get_ylim()[0]),fontsize=14, ha='right', fontstyle='italic')

    if not show_legend:
        ax.legend().remove()

    return kruskal_pval

def plot_umap_pca_grid(X, meta_s, figsize=(34,8), cmap=CPTAC_CMAP, title=None, **kwargs):
    """
    Umap-PCA Grid.

    Args:
        * X: (m feats x n samples)
        * meta_s: pd.Series to colorby
        * figsize
        * cmap: colormap for meta_s
    """
    import matplotlib.gridspec as gridspec

    p_df, pca, feats = get_pca_umap(X, **kwargs)

    fig = plt.figure(constrained_layout=False, figsize=figsize)

    spec = gridspec.GridSpec(nrows=2, ncols=8, wspace=0.3)

    umap_ax = fig.add_subplot(spec[:2,:2])

    pca1_ax = fig.add_subplot(spec[0,2])
    pca2_ax = fig.add_subplot(spec[1,2])
    pca3_ax = fig.add_subplot(spec[0,3])
    pca4_ax = fig.add_subplot(spec[1,3])

    meta_s = meta_s.loc[X.columns]

    _ = plot_umap(p_df, meta_s, s=30, cdict=cmap, ax=umap_ax)
    _ = plot_pca_ax(p_df, pca, cohort_s=meta_s, cohort_colors=cmap, ax=pca1_ax, order=[1,2,3])
    _ = plot_pca_ax(p_df, pca, cohort_s=meta_s, cohort_colors=cmap, ax=pca2_ax, order=[2,3,4])
    _ = plot_pca_ax(p_df, pca, cohort_s=meta_s, cohort_colors=cmap, ax=pca3_ax, order=[3,4,5])
    _ = plot_pca_ax(p_df, pca, cohort_s=meta_s, cohort_colors=cmap, ax=pca4_ax, order=[4,5,6])

    pca1_ax.legend().remove()
    pca2_ax.legend().remove()
    pca3_ax.legend().remove()
    pca4_ax.legend().remove()

    umap_ax.legend(bbox_to_anchor=(-.2, 0.5), loc='center left')

    if title is not None:
        umap_ax.set_title(title, fontsize=20, ha='left', x=0)

    return fig, p_df, pca, feats

# -------------------------
# X
# -------------------------
def plot_cell_fraction(
    cf: pd.DataFrame,
    order: Union[None, pd.Series] = None,
    figsize: tuple = (12,4),
    display_avg: bool = False,
    ax = None,
    vlines: list = None,
    cmap: str = "Set2"
    ):
    """
    Plot cell fraction from deconvolution.
    ------------------------
    Args:
        * cf: pd.DataFrame of deconvolution results
        * order: ordering of x-axis
        * figsize: size of figure
        * display_avg: display average proportion for each cell-type
        * ax: matplotlib axis
        * vlines: list of x-liens
        * cmap: colormap to use for plots

    Returns:
        * axis

    """
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    if order is None:
        cf.sort_values(by=[cf.mean(0).idxmax(axis=0)],ascending=False).plot(
            kind='area',
            stacked=True,
            ax=ax,
            rasterized=True,
            colormap=cmap,
            linewidth=0,
        )
    else:
        missing_samples = list(set(order)-set(cf.index))
        cf_missing = pd.DataFrame(None, index=missing_samples, columns=cf.columns)

        cf = pd.concat((cf,cf_missing))

        cf.loc[order].plot(
            kind='area',
            stacked=True,
            ax=ax,
            rasterized=True,
            colormap=cmap,
            linewidth=0,
        )

    handles, labels = ax.get_legend_handles_labels()
    labels.reverse()
    handles.reverse()

    mean_values = cf.mean(0)[labels]

    if display_avg:
        labels = [x+' - '+"{:.2f}".format(y) for x,y in zip(labels,mean_values)]

    ax.legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.65))

    ax.set_xlim([0,cf.shape[0]-1])
    ax.set_ylim([0,1])
    ax.set_ylabel('Proportion', fontsize=14)
    ax.set_xlabel('Samples', fontsize=14)
    ax.set_xticklabels('')

    if vlines is not None:
        for x in vlines:
            ax.axvline(x, 0, 1, color='black', alpha=0.5)

    return ax

def plot_strip(
    X: pd.DataFrame,
    group_s: pd.Series,
    yax: str,
    figsize: tuple = (10,3),
    s: float = 5,
    order: Union[None, list] = None,
    ylabel: Union[None, str] = None,
    plot_violin: bool = True,
    plot_mean_line: bool = True,
    y_offset: float = 0.5,
    ylim: Union[None, tuple] = None,
    show_pval: Union[str,None] = None,
    ax = None,
    format_fig: bool = False,
    show_legend: bool = False,
    annot = None,
    color = None,
    **kwargs
    ):
    """
    Plot Strip Plot.
    -----------------------
    Args:
        * X: pd.DataFrame (features x samples)
        * group_s: pd.Series with annotations for samples
        * yax: feature to plot
        * figsize: tuple for figure size
        * s: size of dots
        * order: ordering of subtypes
        * ylabel: str
    Returns:
        * figure
    """
    from scipy.stats import kruskal

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)

    categ = list(set(group_s))
    pos = [x+1 for x in range(len(categ))]

    split_by = group_s.name
    _X = X.T.join(group_s)
    _X[yax] = _X[yax].astype(float)

    kruskal_statistic, kruskal_pval = kruskal(*_X.groupby(split_by)[yax].apply(np.array).values)

    if color is None:
        color = sns.color_palette("husl", len(categ))

    sns.stripplot(
        x=split_by,
        y=yax,
        data=_X,
        ax=ax,
        s=s,
        palette=color,
        order=order,
        **kwargs
    )

    if plot_mean_line:
        # Add means
        capprops = dict(linewidth=0, color='white')
        sns.boxplot(
            x=split_by,
            y=yax,
            data=_X,
            ax=ax,
            showfliers=False,
            showbox=False,
            whis=0,
            capprops=capprops,
            order=order
        )

    if plot_violin:
        sns.violinplot(
            x=split_by,
            y=yax,
            data=_X,
            ax=ax,
            alpha=0.2,
            color='white',
            order=order
        )

    if ylim is None:
        ax.set_ylim((ax.get_ylim()[0]-y_offset,ax.get_ylim()[1]+y_offset))
    else:
        ax.set_ylim(ylim)

    # Y label
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    # Add axis
    leg = ax.legend([yax], handlelength=0, handletextpad=0, fancybox=True, fontsize=12)
    for item in leg.legendHandles:
        item.set_visible(False)

    # Add p-value
    if show_pval is not None:
        if show_pval == 'lower left':
            ax.annotate("   Kruskal-Wallis, p = {:.3E}\n".format(kruskal_pval),(ax.get_xlim()[0],ax.get_ylim()[0]),fontsize=8)
        elif show_pval == 'above left':
            ax.annotate("   Kruskal-Wallis, p = {:.3E}\n".format(kruskal_pval),(ax.get_xlim()[0],ax.get_ylim()[1]),fontsize=8)

    if format_fig:
        _format_ax(ax)

    if annot is not None:
        ax.annotate("{}\n".format(annot),(ax.get_xlim()[1], ax.get_ylim()[0]),fontsize=14, ha='right', fontstyle='italic')

    if not show_legend:
        ax.legend().remove()

    return kruskal_pval

def _plot_matrix_heatmap(mat, ax=None, cbar_ax=None, vmin=None, vmax=None, auto_thresh=False, cmap='coolwarm', **kwargs):
    """
    Plot pathway matrix.
    """
    if ax is None:
        fig,ax = plt.subplots(figsize=(8,4))
        cbar_ax = fig.add_axes(
            [-ax.get_position().x1*0.1, ax.get_position().y0, .05, ax.get_position().y1-ax.get_position().y0]
        )
        cbar_ax.yaxis.set_ticks_position('left')

    if auto_thresh:
        vmax = max(mat.max().max(), -mat.min().min())
        vmin = -vmax

    sns.heatmap(
        mat.T,
        ax=ax,
        cbar_ax=cbar_ax,
        vmin=vmin,
        vmax=vmax,
        cbar=False,
        cmap=cmap,
        **kwargs
    )

    ax.hlines(range(mat.shape[1]+1), *ax.get_xlim(), linewidth=0.5)

    ax.set_xticks([])
    ax.set_yticks([x+.5 for x in np.arange(mat.shape[1])])
    ax.set_yticklabels(mat.columns)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.yaxis.tick_right()
    ax.tick_params(length=0)
    ax.set_yticklabels(ax.get_yticklabels(), rotation=0)

    for _, spine in ax.spines.items():
        spine.set_visible(True)

def plot_cluster_heatmap(
    X: pd.DataFrame,
    metric: str = 'euclidean',
    method: str = 'ward',
    n_clusters: int = 10,
    color_thresh_scale: float = 0.3,
    figsize: tuple = (8,8),
    p: int = 30,
    vmax: Union[float, None] = None,
    vmin: Union[float, None] = None,
    cbar_label: str = 'Z',
    cmap: Union[str, None] = None,
    show_h_lines: bool = True,
    hms = None,
    cluster_cdict = None,
    cluster_names = None,
    cohort_m = None,
    cluster_y = True,
    cluster_x = True,
    show_top_dendrogram = True
    ):
    """
    Plot cluster heatmap.
    -----------------------
    Args:
        * X: input matrix (features x samples)
        * metric: distance metric
        * method: method of clustering
        * n_clusters: number of clusters for agglomerative clustering
        * color_thresh_scale: asthetic scale for coloring of dendrogram
        * figsize: figsize
        * p: parameter for dendrogram
        * vmax
        * vmin
        * cbar_label: label for cbar
        * cmap: colormap
        * show_h_lines: separate categories on y axis

    Returns:
        * fig, clust_df: pd.DataFrame of results, feature_order
    """
    from scipy.cluster import hierarchy
    import scipy.cluster.hierarchy as shc
    from sklearn.cluster import AgglomerativeClustering
    from matplotlib.pyplot import cm
    import matplotlib as mpl

    # -------------
    # Heatmap
    # -------------
    fig,ax = plt.subplots(figsize=figsize)
    cbar_ax = fig.add_axes([ax.get_position().x0-.15, ax.get_position().y1+0.05, .025, .1])

    # Compute initial linkage to grab ordering
    d_linkage = shc.linkage(X.T, metric=metric, method=method)
    dres = shc.dendrogram(d_linkage, p=p, no_plot=True)

    if cluster_x:
        dgram_idx = list(map(int, dres['ivl']))
    else:
        dgram_idx = np.arange(X.shape[1])

    # Compute initial linkage to grab ordering
    if cluster_y:
        d_linkage_y = shc.linkage(X, metric=metric, method=method)
        dres_y = shc.dendrogram(d_linkage_y, p=p, no_plot=True, optimal_ordering=True)
        dgram_idx_y = list(map(int, dres_y['ivl']))
    else:
        dgram_idx_y = np.arange(X.shape[0])

    # Create heatmap
    if vmax is None:
        cbar_top_lim = np.max(X.values)
    else:
        cbar_top_lim = vmax

    if vmin is None:
        cbar_bottom_lim = 0
    else:
        cbar_bottom_lim = vmin

    # Create heatmap
    sns.heatmap(
        X.iloc[dgram_idx_y[::-1], dgram_idx].values,
        ax=ax,
        cbar_ax=cbar_ax,
        cbar_kws = {'ticks':[cbar_bottom_lim, cbar_top_lim]},
        rasterized=True,
        vmax=vmax,
        vmin=vmin,
        cmap=cmap
    )

    cbar_ax.set_ylabel(cbar_label, fontsize=10,rotation=90)
    ax.set_xticks([])

    x0 = ax.get_position().x0
    x1 = ax.get_position().x1
    y0 = ax.get_position().y0
    y1 = ax.get_position().y1

    buf = y1*0.015

    # -------------
    # Clustering
    # -------------
    cluster = AgglomerativeClustering(
        n_clusters=n_clusters,
        affinity=metric,
        linkage=method
    )

    clusters = cluster.fit_predict(X.iloc[dgram_idx_y[::-1], dgram_idx].T)
    cluster_color_list, colors = sa.pl.series_to_colors(pd.Series(clusters), cdict=cluster_cdict)

    # -------------
    # Top Dendrogram
    # -------------
    if show_top_dendrogram:
        cmap = cm.rainbow(np.linspace(0, 1, 10))
        hierarchy.set_link_color_palette([mpl.colors.rgb2hex(rgb[:3]) for rgb in cmap])

        dax = fig.add_axes([x0, y1+2*buf, x1-x0, 0.15])

        dres = shc.dendrogram(
            d_linkage,
            p=p,
            ax=dax,
            above_threshold_color="grey",
            color_threshold=color_thresh_scale*max(d_linkage[:,2])
        )

        dax.set_xticks([])
        dax.set_yticks([])
        [dax.spines[x].set_visible(False) for x in ['top','right','bottom','left']]

    # -------------
    # Clusters Axis
    # -------------
    tax = fig.add_axes([x0, y1+.5*buf, x1-x0, buf])
    mat, cmap = sa.pl.color_list_to_matrix_and_cmap(cluster_color_list, order_dict={v: k for k, v in colors.items()})
    tax_cax = fig.add_axes([x1+.5*buf, y1+.5*buf, buf, n_clusters*buf*1.5])

    sns.heatmap(
        mat,
        cmap=cmap,
        ax=tax,
        xticklabels=False,
        yticklabels=False,
        rasterized=False,
        cbar_ax=tax_cax,
        cbar_kws=dict(ticks=range(n_clusters)),
        norm=mpl.colors.BoundaryNorm(np.array(range(n_clusters+1))-0.5, n_clusters)
    )

    if cluster_names is not None:
        tax_cax.set_yticklabels([cluster_names[int(x.get_text())] for x in tax_cax.get_yticklabels()])

    [tax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]
    [tax_cax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

    # -------------
    # Left Dendrogram
    # -------------
    if cluster_y:
        cmap = cm.rainbow(np.linspace(0, 1, 10))
        lax = fig.add_axes([x0-buf*.25-0.07, y0, 0.07, y1-y0])

        dres_lax = shc.dendrogram(
            d_linkage_y,
            p=p,
            ax=lax,
            above_threshold_color="grey",
            orientation="left",
            color_threshold=0
        )

        lax.set_xticks([])
        lax.set_yticks([])
        [lax.spines[x].set_visible(False) for x in ['top','right','bottom','left']]

    # -------------
    # Feature Labels & final Clust
    # -------------
    y_labels = X.iloc[dgram_idx_y[::-1], dgram_idx].index
    ax.yaxis.tick_right()
    ax.tick_params(axis=u'both', which=u'both',length=0)
    ax.set_yticklabels(y_labels, rotation=0)

    clust_df = pd.DataFrame(clusters, index=X.columns[dgram_idx]).rename(columns={0:'clusters'})
    clust_df['c'] = cluster_color_list.values

    if show_h_lines:
        ax.hlines(range(X.shape[1]), *ax.get_xlim(), linewidth=0.5)
    [ax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

    # Add rows underneath
    if hms is not None:
        n_size1 = 5
        bax = fig.add_axes([x0, y0-2.5*buf, x1-x0, 2*buf])

        _cohort_s,_ = sa.pl.series_to_colors(cohort_m.iloc[dgram_idx], cdict=CPTAC_CMAP)
        _cohort_m,_cohort_cmap = sa.pl.color_list_to_matrix_and_cmap(_cohort_s)
        sns.heatmap(_cohort_m,cmap=_cohort_cmap,ax=bax,
            xticklabels=False,yticklabels=False,rasterized=True,cbar=False
        )

        bax.set_yticks([0.5])
        bax.set_yticklabels(['Cohort'])
        bax.yaxis.tick_right()
        bax.tick_params(length=0)

        [bax.spines[x].set_visible(True) for x in ['top','right','bottom','left']]

        bax = fig.add_axes([x0, y0-n_size1*buf, x1-x0, (n_size1-3)*buf])
        _plot_matrix_heatmap(hms[0].iloc[dgram_idx], ax=bax, auto_thresh=True)

        n_size = 12
        bax = fig.add_axes([x0, y0-n_size1*buf-n_size*buf, x1-x0, (n_size-1)*buf])
        _plot_matrix_heatmap(hms[1].iloc[dgram_idx], ax=bax, cmap="Reds")

    return fig, clust_df, np.array(X.index[dgram_idx_y][::-1])


# -------------------------
# UMAP
# -------------------------
def get_umap(H):
    """
    Compute UMAP. Add to H-matrix.
    """
    import umap

    assert 'umap1' not in list(H), 'UMAP already computed.'
    fit = umap.UMAP()
    H_umap = fit.fit_transform(H.iloc[:,:-3])
    H['umap1'] = H_umap[:,0]
    H['umap2'] = H_umap[:,1]

def plot_umap(H, c, cdict=None, on_plot=False, figsize=(8,8), s=10, alpha=0.8, ax=None, **kwargs):
    """
    Plots single UMAP of H-matrix.
    Args
        H: H-matrix result from NMF (pandas dataframe - n_samples x n_signatures)
        c: pd.Series
        cdict: color dict
        on_plot: whether or not to plot on the plot
        figsize: figure size
        s: size of points
        alpha: alpha
    Returns
        Matplotlib Figure
    """
    import matplotlib.patheffects as path_effects

    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    if isinstance(cdict,dict):
        subtypes_cdict = cdict
    else:
        if cdict is None:
            cdict='husl'
        subtypes = list(set(c))
        subtypes_cdict = {subtypes[i]:x for i,x, in enumerate(sns.color_palette(cdict, len(set(c))))}

    sc = ax.scatter(
        H['umap1'].values,
        H['umap2'].values,
        c=c[H.index].apply(lambda x: np.array(subtypes_cdict[x])),
        alpha=alpha,
        s=s,
        rasterized=True,
        edgecolor='black',
        linewidth=0.1
    )

    ax.set_xlabel(r'UMAP1 $\to$', fontsize=18)
    ax.xaxis.set_label_coords(0.12,-0.025)
    ax.set_ylabel(r'UMAP2 $\to$', fontsize=18)
    ax.yaxis.set_label_coords(-0.02,0.12)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlim(*ax.get_xlim())
    ax.set_ylim(*ax.get_ylim())

    if on_plot:
        # Annotate each cluster
        for lab in set(c):
            txt = ax.annotate(lab, \
                   H.loc[c[c==lab].index,['umap1','umap2']].mean(), \
                   horizontalalignment='center', \
                   verticalalignment='center', \
                   size=16, \
                   color='white' \
            )

            txt.set_path_effects([path_effects.Stroke(linewidth=2, foreground='black'), path_effects.Normal()])
    else:
        for k,v in subtypes_cdict.items():
            m = ax.scatter(-1e2, -1e2, alpha=.8, c=v, label=k)

        ax.legend(**kwargs)

def plot_umap_sigs(H, H_rep, ncols=4, s=3, alpha=0.7, width_m=5, height_m=4):
    """
    Plot Signatures.
    """
    ids = list(set(H['max_id']))
    nrows = int(np.ceil((len(ids))/ncols))

    fig,axes = plt.subplots(nrows, ncols, figsize=(ncols*width_m, nrows*height_m))

    ax_idx = 0
    for row in range(nrows):
        for col in range(ncols):
            try:
                im = axes[row,col].scatter(
                        H_rep[:,0], \
                        H_rep[:,1], \
                        c=H['S{}'.format(ids[ax_idx])], \
                        s=s, \
                        alpha=alpha, \
                        cmap='Reds', \
                        rasterized=True \
                )

                axes[row,col].set_title('S{}'.format(ids[ax_idx]))
                axes[row,col].set_xticks([])
                axes[row,col].set_yticks([])

                if col==0 and row==(nrows-1):
                    axes[row,col].set_xlabel(r'UMAP1 $\to$', fontsize=10)
                    axes[row,col].xaxis.set_label_coords(0.12,-0.025)
                    axes[row,col].set_ylabel(r'UMAP2 $\to$', fontsize=10)
                    axes[row,col].yaxis.set_label_coords(-0.02,0.12)

                plt.colorbar(im,ax=axes[row,col])

            except:
                axes[row,col].axis('off')

            ax_idx+=1

    return fig

def plot_dendro_umap(X, group_s, ax=None, figsize=(5,5)):
    """
    Plot dendro umap.
    """
    if ax is None:
        fig,ax = plt.subplots(figsize=figsize)

    _ = plot_umap(
        X,
        group_s,
        s=30,
        cdict={True:'lightpink', False:'lightblue', np.nan:'white'},
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        frameon=False,
        ax=ax
    )
    ax.legend().remove()
    ax.set_title("Dendro Group {}".format(group_s.name))

def plot_dendro_umap_grid(X, groups_df, ncols=4, s=3, alpha=0.7, width_m=4, height_m=4):
    """
    Plot Signatures.
    """
    ids = list(groups_df.columns)
    nrows = int(np.ceil((len(ids))/ncols))

    fig,axes = plt.subplots(nrows, ncols, figsize=(ncols*width_m, nrows*height_m))

    ax_idx = 0
    for row in range(nrows):
        for col in range(ncols):
            try:
                plot_dendro_umap(X, groups_df[ids[ax_idx]], ax=axes[row,col])

                if col==0 and row==(nrows-1):
                    axes[row,col].set_xlabel(r'UMAP1 $\to$', fontsize=10)
                    axes[row,col].xaxis.set_label_coords(0.12,-0.025)
                    axes[row,col].set_ylabel(r'UMAP2 $\to$', fontsize=10)
                    axes[row,col].yaxis.set_label_coords(-0.02,0.12)
                else:
                    axes[row,col].set_xlabel("")
                    axes[row,col].set_ylabel("")
            except:
                axes[row,col].axis('off')

            ax_idx+=1

    return fig

# -------------------------
# Correlations
# -------------------------
def corr_clustermap(
    corr_df,
    comp,
    a_geneset,
    b_geneset,
    filter_pval=False,
    vmin=-0.7,
    vmax=0.7,
    a_isPTM=False,
    b_isPTM=False,
    rotate=False,
    figsize=None,
    cmap='coolwarm',
    **kwargs
    ):
    """
    Create correlation heatmap.
    -------------------------
    Args:
        * corr_df: dataframe of correlations
        * comp: comparison column (ex. rna_prot)
        * a_geneset: geneset to subset
        * b_geneset: geneset to subset
        * filter_pval: mask p-values
    """
    corr_sub_df = corr_df[
        (corr_df['comp']==comp)
        & (corr_df['a_geneset']==a_geneset)
        & (corr_df['b_geneset']==b_geneset)
    ]

    # Rename PTM sites to be GeneName_AminoAcid
    if a_isPTM:
        corr_sub_df['a_gene'] = corr_sub_df.apply(lambda x: x['a_gene'] + '|' + x['a'].split('_')[2],1)
    if b_isPTM:
        corr_sub_df['b_gene'] = corr_sub_df.apply(lambda x: x['b_gene'] + '|' + x['b'].split('_')[2],1)

    corr_mat = corr_sub_df.pivot(index=['a_gene'], values='rho', columns=['b_gene']).astype(float)

    if filter_pval:
        corr_f_df = corr_sub_df[corr_sub_df['pval']<.1]
        corr_f_df = corr_f_df.pivot(index=['a_gene'], values='rho', columns=['b_gene']).astype(float)
        mask = corr_f_df.isna()
        corr_mat = corr_mat.loc[mask.index, mask.columns]
    else:
        mask = None

    corr_sub_df = corr_sub_df.pivot(index=['a_gene'], values='rho', columns=['b_gene']).astype(float)

    if rotate:
        corr_sub_df = corr_mat.T
        if isinstance(mask, pd.DataFrame):
            mask = mask.T

    if figsize is None:
        figsize=(
            int(corr_sub_df.shape[0]/5),
            int(corr_sub_df.shape[1]/5)
        )

    g = sns.clustermap(
        corr_sub_df,
        yticklabels=True,
        xticklabels=True,
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
        mask=mask,
        figsize=figsize,
        **kwargs
    )

    if rotate:
        g.ax_heatmap.set_xlabel(comp.split("_")[0].upper() + " | " + a_geneset, fontsize=24)
        g.ax_heatmap.set_ylabel(comp.split("_")[1].upper() + " | " + b_geneset, fontsize=24)
    else:
        g.ax_heatmap.set_xlabel(comp.split("_")[1].upper() + " | " + b_geneset, fontsize=24)
        g.ax_heatmap.set_ylabel(comp.split("_")[0].upper() + " | " + a_geneset, fontsize=24)
    return g

# --------------------------------------------------------------------------------------------
# Author: Francois Auget
# --------------------------------------------------------------------------------------------
def setup_figure(aw=4.5, ah=3, xspace=[0.75,0.25], yspace=[0.75,0.25],
                 colorbar=False, ds=0.15, cw=0.15, ct=0, ch=None):
    dl, dr = xspace
    db, dt = yspace
    fw = dl + aw + dr
    fh = db + ah + dt
    fig = plt.figure(facecolor=(1,1,1), figsize=(fw,fh))
    ax = fig.add_axes([dl/fw, db/fh, aw/fw, ah/fh])
    if not colorbar:
        return ax
    else:
        if ch is None:
            ch = ah/2
        cax = fig.add_axes([(dl+aw+ds)/fw, (db+ah-ch-ct)/fh, cw/fw, ch/fh])
        return ax, cax

def hexdensity(x, y, bounds=None, bins='log', scale='log',
               cmap=None, vmin=None, vmax=None, ax=None, cax=None,
               unit='TPM', entity='genes',
               gridsize=175, fontsize=12, show_corr=True, clip_on=True, rasterized=False):
    """Wrapper for hexbin"""
    import scipy
    import matplotlib.ticker as ticker
    if ax is None: # setup new axes
        ax, cax = setup_figure(2, 2, xspace=[0.75, 1], yspace=[0.75, 0.5], colorbar=True, ch=1, cw=0.12)
        ax.margins(0.01)

    if cmap is None:
        cmap = copy.copy(plt.cm.RdYlBu_r)
        cmap.set_bad('w', 1.)

    rho = scipy.stats.spearmanr(x, y)[0]
    x = x.copy()
    y = y.copy()
    nanidx = (x == 0) | (y == 0)
    x[nanidx] = np.NaN
    y[nanidx] = np.NaN

    h = ax.hexbin(x, y, bins=bins, xscale=scale, yscale=scale, linewidths=0.1,
                  gridsize=gridsize, cmap=cmap, vmin=vmin, vmax=vmax, mincnt=1, zorder=1,
                  clip_on=clip_on, rasterized=rasterized)

    if bounds is None:
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        bounds = [np.minimum(xlim[0], ylim[0]), np.maximum(xlim[1], ylim[1])]
    elif len(bounds) == 2:
        ax.set_xlim(bounds)
        ax.set_ylim(bounds)
    else:
        ax.set_xlim(bounds[:2])
        ax.set_ylim(bounds[2:])
    ax.spines['left'].set_position(('outward', 6))
    ax.spines['bottom'].set_position(('outward', 6))

    if show_corr:
        t = ax.text(.95, 0, r'$\rho$ = {:.2f}'.format(rho), transform=ax.transAxes,
                    ha='right', va='bottom', fontsize=fontsize, zorder=2)
        t.set_bbox(dict(facecolor='w', alpha=0.5, edgecolor='none', boxstyle="round,pad=0.1"))

    #hc = plt.colorbar(h, cax=cax, orientation='vertical', ticks=ticker.LogLocator(numticks=4))
    #hc.set_label('log$\mathregular{_{10}}$('+entity+')', fontsize=fontsize)

    if isinstance(x, pd.Series):
        ax.set_xlabel(f'{x.name} ({unit})', fontsize=fontsize)
    if isinstance(y, pd.Series):
        ax.set_ylabel(f'{y.name} ({unit})', fontsize=fontsize)

    return ax, cax

def qqplot(pval, pval_null=None, title='', labels=None, ax=None, c=None, s=16):
    """QQ-plot"""
    if labels is None:
        labels = ['', '']

    n = len(pval)
    x = -np.log10(np.arange(1,n+1)/(n+1))

    if ax is None:
        ax = setup_figure(4,4)

    ax.margins(x=0.02, y=0.05)
    args = {'s':s, 'edgecolor':'none', 'clip_on':False, 'alpha':1, 'rasterized':True}
    log_pval_sorted = -np.log10(np.sort(pval))

    ax.scatter(
        x,
        log_pval_sorted,
        c=c,
        zorder=30,
        label=labels[0],
        **args
    )

    if pval_null is not None:
        assert len(pval)==len(pval_null)
        log_pval_sorted = -np.log10(np.sort(pval_null))
        ax.scatter(
            x,
            log_pval_sorted,
            c=[[0.5]*3],
            zorder=20,
            label=labels[1],
            **args
        )

    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, min_n_ticks=5, nbins=4))
    ax.set_xlabel('Expected -log$\mathregular{_{10}}$(p-value)', fontsize=14)
    ax.set_ylabel('Observed -log$\mathregular{_{10}}$(p-value)', fontsize=14)

    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    ax.set_xlim([0, xlim[1]])
    ax.set_ylim([0, ylim[1]])
    ci = 0.95
    xi = np.arange(1, n+1)
    clower = -np.log10(stats.beta.ppf((1-ci)/2, xi, xi[::-1]))
    cupper = -np.log10(stats.beta.ppf((1+ci)/2, xi, xi[::-1]))
    ax.fill_between(x, cupper, clower, color=[[0.8]*3], clip_on=True, rasterized=True)
    ax.plot([x[0], x[-1]], [x[0], x[-1]], '--', lw=1, color=[0.2]*3, zorder=50, clip_on=True, rasterized=True)
    #ax.spines['left'].set_position(('outward', 6))
    #ax.spines['bottom'].set_position(('outward', 6))
    ax.set_title('{}'.format(title), fontsize=12)
    if labels[0] != '':
        ax.legend(loc='upper left', fontsize=10, handlelength=0.5, handletextpad=0.33)
# --------------------------------------------------------------------------------------------
