import pandas as pd
import numpy as np
import os
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import sklearn.decomposition
import subprocess
import sys
from cmapPy.pandasGEXpress.parse_gct import parse
from subprocess import PIPE
from collections import defaultdict
from scipy import stats
from qtl.norm import deseq2_size_factors

def tpm_loader(tpm, counts, samples=None, impute_zero_to_min=False):
    """
    Bulk load dataset.
    """
    # Load data
    if tpm.endswith(".gct"):
        tpm = pd.read_csv(tpm, sep='\t', skiprows=2, index_col=0)
    else:
        tpm = pd.read_parquet(tpm)

    if counts.endswith(".gct"):
        counts = pd.read_csv(counts, sep='\t', skiprows=2, index_col=0)
    else:
        counts = pd.read_parquet(counts)

    gene_name = tpm.loc[:,['Description']]
    tpm = tpm.iloc[:,1:]

    if samples is not None:
        tpm = tpm.loc[:,samples]

    # Filter counts
    tpm = tpm[(np.sum(tpm >= 0.1, 1) > tpm.shape[1]*0.2) & (np.sum(counts.iloc[:,1:] >= 6, 1) > tpm.shape[1]*0.2)]

    # Impute
    if impute_zero_to_min:
        tpm_imp = tpm.copy()
        for x in tpm_imp.values:
            m = x == 0
            if any(m):
                x[m] = min(x[~m]) / 2
        return tpm_imp, np.log2(tpm_imp / deseq2_size_factors(tpm_imp)), counts, gene_name
    else:
        return tpm, np.log2(1+tpm / deseq2_size_factors(tpm)), counts, gene_name

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

# Plot PCA Genes
def plot_pca_genes(
    P_df,
    pca,
    genes,
    tpm,
    order=[1,2],
    vmin=None,
    vmax=None,
    ncols=2,
    s=3,
    alpha=0.7,
    width_m=4,
    height_m=3,
    cmap=plt.cm.Spectral_r
    ):
    """
    Plots a set of PC genes.
    """
    nrows = int(np.ceil((len(genes)/ncols)))
    fig,axes = plt.subplots(nrows, ncols, figsize=(ncols*width_m, nrows*height_m))

    ax_idx = 0
    for row in range(nrows):
        for col in range(ncols):
            try:
                im = axes[row,col].scatter(
                    P_df[order[1]-1],
                    P_df[order[0]-1],
                    c=tpm.loc[genes[ax_idx], P_df.index],
                    cmap=cmap,
                    vmin=vmin,
                    vmax=vmax,
                    alpha=alpha,
                    s=s,
                    rasterized=True
                )

                axes[row,col].set_title('{}'.format(genes[ax_idx]), fontsize=14)
                axes[row,col].set_xticks([])
                axes[row,col].set_yticks([])
                plt.colorbar(im,ax=axes[row,col])

            except:
                axes[row,col].axis('off')

            ax_idx+=1

    # Set Axis Labels
    fig.text(0.5, 0.08, 'PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), ha='center', va='center', fontsize=14)
    fig.text(0.08, 0.5, 'PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), ha='center', va='center', rotation='vertical', fontsize=14)

    return fig

def plot_pca_ax(
    P_df,
    pca,
    ax=None,
    c=None,
    cohort_s=None,
    cohort_colors=None,
    cohort_args=None,
    order=[1,2,3],
    outliers=None,
    title='',
    vmin=None,
    vmax=None,
    alpha=1,
    lw=0,
    s=30,
    cmap=plt.cm.Spectral_r,
    cticks=None,
    cticklabels=None,
    clabel='',
    show_legend=True
    ):
    """
    PCA Plot by axis.
    -------------------
    cohort_s: Series encoding cohorts
    cohort_colors: dict
    Modes:
    """
    if cohort_s is not None:
        cohorts = cohort_s.unique()
        nc = len(cohorts)
        if cohort_colors is None and cohort_args is None:
            cohort_colors = {i:j for i,j in zip(cohorts, sns.husl_palette(nc, s=1, l=0.6))}
        if cohort_args is None:
            cohort_args = {}
            for k in np.unique(cohort_s):
                cohort_args[k] = {'color': cohort_colors[k], 'marker':'o', 'edgecolor':'none', 's':s}

    if ax is None:
        fig,ax = plt.subplots(figsize=(6,6))

    if cohort_s is None:
        sa = ax.scatter(P_df[order[1]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s, rasterized=False)
    else:
        for k in np.unique(cohort_s):
            i = cohort_s[cohort_s==k].index
            ax.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, rasterized=True, **cohort_args[k])

    format_plot(ax, fontsize=10)
    ax.set_xlabel('PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), fontsize=12)
    ax.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if outliers is not None:
        ax.scatter(P_df.loc[outliers, order[1]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None, rasterized=True)

    ax.set_title(title, fontsize=12)

    if cohort_s is not None and show_legend:
        leg = ax.legend(loc=0, fontsize=9, scatterpoints=1, handletextpad=0.1, framealpha=1, labelspacing=0.35)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    if cohort_s is None and c is not None and len(c)==P_df.shape[0]:
        x1 = ax.get_position().x1
        y1 = ax.get_position().y1

        fig = plt.gcf()
        cax = fig.add_axes(np.array([x1, y1*5/5.5, 0.15/5.5, 1/5.5]))
        hc = plt.colorbar(sa, cax=cax, orientation='vertical')

        if cticks is not None:
            hc.set_ticks(cticks)
        if cticklabels is not None:
            hc.ax.tick_params(labelsize=9)
            cax.set_xticklabels(cticklabels, fontsize=10)

        hc.locator = ticker.MaxNLocator(integer=True, min_n_ticks=2, nbins=5)
        hc.update_ticks()

        cax.set_ylabel(clabel, rotation=0, ha='right', va='center', fontsize=12)

    return ax

def plot_pca_gridplot(X, meta, norm_sites, title=None, c='cohort', figsize=(16,12), cohort_colors=None):
    """
    PCA Grid Plot.
    ---------------
    Built on Francois's PCA plotting code.
    Plots a grid of PCA's for each PTM colored by cohort.

    Args:
        * X: dict of numpy matrices
        * meta: dict of metadata entries
        * norm_sites: "stable" sites across all samples
        * title: title of plot
        * c: variable to color scatterplots by (must be in meta[key].columns)
        * figsize: tuple of figure size

    Returns:
        * matplotlib Figure
    """
    fig,axes = plt.subplots(3,4, figsize=figsize)

    for f_idx, feat in enumerate(X.keys()):
        # All Sites
        P_df, pca, prots = get_pcs(X[feat].T, normalize=False, return_genes=True)

        # PC 1,2,3
        _ = plot_pca_ax(P_df, pca, ax=axes[f_idx,0], order=[1,2], cohort_s=meta[feat].loc[X[feat].index][c].sort_values(), cohort_colors=cohort_colors)
        _ = plot_pca_ax(P_df, pca, ax=axes[f_idx,1], order=[1,3], cohort_s=meta[feat].loc[X[feat].index][c].sort_values(), cohort_colors=cohort_colors)

        axes[f_idx,1].set_ylabel("")
        axes[f_idx,1].set_yticks([])
        axes[f_idx,1].spines['left'].set_linewidth(0)

        # Stable Sites
        P_df, pca, prots = get_pcs(X[feat].T.loc[norm_sites[feat]], normalize=False, return_genes=True)

        # PC 1,2,3
        _ = plot_pca_ax(P_df, pca, ax=axes[f_idx,2], cohort_s=meta[feat].loc[X[feat].index][c].sort_values(), cohort_colors=cohort_colors)
        _ = plot_pca_ax(P_df, pca, ax=axes[f_idx,3], cohort_s=meta[feat].loc[X[feat].index][c].sort_values(), cohort_colors=cohort_colors)

        axes[f_idx,3].set_ylabel("")
        axes[f_idx,3].set_yticks([])
        axes[f_idx,3].spines['left'].set_linewidth(0)

        [axes[f_idx,x].get_legend().remove() for x in (0,1,2)]

        pad=5
        axes[f_idx,0].annotate(
            feat.capitalize(),
            xy=(1.3, 1.8),
            xytext=(-axes[f_idx,0].yaxis.labelpad - pad, 0),
            xycoords=axes[f_idx,0].yaxis.label,
            textcoords='offset points',
            fontsize=18,
            weight='bold'
        )

    axes[0,0].set_title(title, fontsize=24, x=0.2, y=1.2)
    axes[0,1].set_title("All Sites", fontsize=18, y=1.1, ha='left')
    axes[0,2].set_title("Stable Sites", fontsize=18, y=1.1, ha='right')

    axes[1,3].get_legend().remove()
    axes[2,3].get_legend().remove()

    plt.tight_layout()

    return fig

#------------------------------------------------------------------
# From Francois
#------------------------------------------------------------------
def remove_covariates(df, C, center=False, fail_colinear=False):
    """
    Residualizes rows of M relative to columns of C
    """

    # transform input
    if isinstance(df, pd.DataFrame) or isinstance(df, pd.Series):
        M = df.values
    else:
        M = df

    isvector = False
    if isinstance(M, list) or (hasattr(M, 'shape') and len(M.shape)==1):
        M = np.array(M).reshape(1,-1)
        isvector = True

    if isinstance(C, list) or (hasattr(C, 'shape') and len(C.shape)==1):
        C = np.array(C).reshape(-1,1)

    Q = orthogonalize_covariates(C, fail_colinear=fail_colinear)

    # residualize M relative to C
    M0 = (M.T - np.mean(M,axis=1)).T
    if center:
        M0 = M0 - np.dot(np.dot(M0, Q), Q.T)
    else:
        M0 = M - np.dot(np.dot(M0, Q), Q.T)  # retain original mean

    if isvector:
        M0 = M0[0]

    if isinstance(df, pd.DataFrame):
        M0 = pd.DataFrame(M0, index=df.index, columns=df.columns)
    elif isinstance(df, pd.Series):
        M0 = pd.Series(M0, index=df.index, name=df.name)

    return M0

def orthogonalize_covariates(C, fail_colinear=True):
    """
    C: covariates (columns)
    """
    # center and orthogonalize
    Q,R = np.linalg.qr(C-np.mean(C,axis=0))

    # check for colinearity
    colinear_ix = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
    if np.any(colinear_ix):
        if fail_colinear:
        # if np.min(np.abs(np.diag(R))) < np.finfo(np.float64).eps * C.shape[1]:
            raise ValueError("Colinear or zero covariates detected")
        else:  # drop colinear covariates
            print('  * Colinear covariates detected. {} covariates dropped.'.format(np.sum(colinear_ix)))
            Q = Q[:, ~colinear_ix]

    return Q

def normalize_counts(gct_df, C=None):
    gct_norm_df = gct_df / deseq2_size_factors(gct_df)
    gct_norm_df = np.log10(1+gct_norm_df)

    # threshold low expressed genes
    mask = np.mean(gct_norm_df > 1, axis=1) > 0.1  # >=10 counts in >10% of samples
    gct_norm_df = gct_norm_df[mask]

    if C is not None:
        gct_norm_df = remove_covariates(gct_norm_df, C, center=False)

    # gct_norm_std_df = center_normalize(gct_norm_df)
    gct_norm_std_df = gct_norm_df - gct_norm_df.mean(axis=0)
    gct_norm_std_df = gct_norm_std_df / np.sqrt(gct_norm_std_df.pow(2).sum(axis=0))

    return gct_norm_std_df

def get_pcs(gct_df, normalize=True, C=None, n_components=5, return_genes=False):
    """
    Scale input GCT, threshold, normalize and calculate PCs
    """
    if normalize:
        gct_norm_std_df = normalize_counts(gct_df, C=C)
    else:
        gct_norm_std_df = gct_df

    pca = sklearn.decomposition.PCA(n_components=n_components)
    pca.fit(gct_norm_std_df.T)
    P = pca.transform(gct_norm_std_df.T)
    P_df = pd.DataFrame(P, index=gct_norm_std_df.columns)

    if return_genes:
        return P_df, pca, gct_norm_std_df.index.values
    else:
        return P_df, pca

def plot_pca(P_df, pca, c=None, cohort_s=None, cohort_colors=None, cohort_args=None, order=[1,2,3], outliers=None, title='',
    vmin=None, vmax=None, alpha=1, lw=0, s=30, cmap=plt.cm.Spectral_r, cticks=None, cticklabels=None, clabel='',
    show_legend=True, show_ax2=True):
    """
    cohort_s: Series encoding cohorts
    cohort_colors: dict
    Modes:
    """
    if cohort_s is not None:
        cohorts = cohort_s.unique()
        nc = len(cohorts)
        if cohort_colors is None and cohort_args is None:
            # cohort_colors = {i:j for i,j in zip(cohorts, cm.get_cmap(cmap, nc)(np.arange(nc)))}
            cohort_colors = {i:j for i,j in zip(cohorts, sns.husl_palette(nc, s=1, l=0.6))}
        if cohort_args is None:
            cohort_args = {}
            for k in np.unique(cohort_s):
                cohort_args[k] = {'color': cohort_colors[k], 'marker':'o', 'edgecolor':'none', 's':s}

    if show_ax2:
        fig = plt.figure(facecolor=(1,1,1), figsize=(10.5,5.5))
        ax1 = fig.add_axes(np.array([1/10.5, 0.75/5.5, 4/10.5, 4/5.5]))
    else:
        fig = plt.figure(facecolor=(1,1,1), figsize=(5.5,5.5))
        ax1 = fig.add_axes(np.array([1/5.5, 0.75/5.5, 4/5.5, 4/5.5]))
    if cohort_s is None:  # c[P_df.index]
        sa = ax1.scatter(P_df[order[1]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
    else:
        for k in np.unique(cohort_s):
        # for k in cohort_s.unique():
            i = cohort_s[cohort_s==k].index
            ax1.scatter(P_df.loc[i,order[1]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
    format_plot(ax1, fontsize=10)
    ax1.set_xlabel('PC {0} ({1:.2f}%)'.format(order[1], pca.explained_variance_ratio_[order[1]-1]*100), fontsize=12)
    ax1.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if show_ax2:
        ax2 = fig.add_axes(np.array([6/10.5, 0.75/5.5, 4/10.5, 4/5.5]))
        if cohort_s is None:
            ax2.scatter(P_df[order[2]-1], P_df[order[0]-1], c=c, cmap=cmap, vmin=vmin, vmax=vmax, lw=lw, alpha=alpha, s=s)
        else:
            for k in np.unique(cohort_s):
                i = cohort_s[cohort_s==k].index
                ax2.scatter(P_df.loc[i,order[2]-1], P_df.loc[i,order[0]-1], alpha=alpha, label=k, **cohort_args[k])
            # ax2.legend(loc=3, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))

        format_plot(ax2, fontsize=10)
        ax2.set_xlabel('PC {0} ({1:.2f}%)'.format(order[2], pca.explained_variance_ratio_[order[2]-1]*100), fontsize=12)
        ax2.set_ylabel('PC {0} ({1:.2f}%)'.format(order[0], pca.explained_variance_ratio_[order[0]-1]*100), fontsize=12)

    if outliers is not None:
        ax1.scatter(P_df.loc[outliers, order[1]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)
        if show_ax2:
            ax2.scatter(P_df.loc[outliers, order[2]-1], P_df.loc[outliers, order[0]-1], c='none', edgecolors='r', marker='s', lw=1, alpha=1, s=50, label=None)

    fig.suptitle(title, fontsize=12)

    if cohort_s is not None and show_legend:
        # ax2.legend(loc=0, fontsize=10, scatterpoints=1, handletextpad=0.1, framealpha=0.5, bbox_to_anchor=(-0.5,-0.1))
        leg = ax1.legend(loc=0, fontsize=9, scatterpoints=1, handletextpad=0.1, framealpha=1, labelspacing=0.35)
        for lh in leg.legendHandles:
            lh.set_alpha(1)

    # if cohort_s is None and c is not None and not isinstance(c, list) and not isinstance(c, str):
    if cohort_s is None and c is not None and len(c)==P_df.shape[0]:
        if show_ax2:
            cax = fig.add_axes(np.array([3.5/10.5, 5/5.5, 1.5/10.5, 0.15/5.5]))
        else:
            cax = fig.add_axes(np.array([3.5/5.5, 5/5.5, 1.5/5.5, 0.15/5.5]))
        # cax = fig.add_axes(np.array([3.5/10.5, 4.85/5.5, 1.5/10.5, 0.15/5.5]))
        hc = plt.colorbar(sa, cax=cax, orientation='horizontal')
        if cticks is not None:
            hc.set_ticks(cticks)
        if cticklabels is not None:
            # hc.set_ticks([0,0.5,1])
            hc.ax.tick_params(labelsize=9)
            # cax.invert_xaxis()
            cax.set_xticklabels(cticklabels, fontsize=10)

        hc.locator = ticker.MaxNLocator(integer=True, min_n_ticks=2, nbins=5)
        hc.update_ticks()

        cax.set_ylabel(clabel, rotation=0, ha='right', va='center', fontsize=12)
    return fig

def format_plot(ax, tick_direction='out', tick_length=4, hide=['top', 'right'], hide_spines=True, lw=1, fontsize=9):

    for i in ['left', 'bottom', 'right', 'top']:
        ax.spines[i].set_linewidth(lw)

    # ax.axis["left"].major_ticklabels.set_ha("left")
    ax.tick_params(axis='both', which='both', direction=tick_direction, labelsize=fontsize)

    # set tick positions
    if 'top' in hide and 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('none')
    elif 'top' in hide:
        ax.get_xaxis().set_ticks_position('bottom')
    elif 'bottom' in hide:
        ax.get_xaxis().set_ticks_position('top')
    else:
        ax.get_xaxis().set_ticks_position('both')

    if 'left' in hide and 'right' in hide:
        ax.get_yaxis().set_ticks_position('none')
    elif 'left' in hide:
        ax.get_yaxis().set_ticks_position('right')
    elif 'right' in hide:
        ax.get_yaxis().set_ticks_position('left')
    else:
        ax.get_yaxis().set_ticks_position('both')

    if hide_spines:
        for i in hide:
            ax.spines[i].set_visible(False)


    # adjust tick size
    for line in ax.xaxis.get_ticklines() + ax.yaxis.get_ticklines():
    #for line in ax.get_xticklines() + ax.get_yticklines():
        line.set_markersize(tick_length) # tick length
        line.set_markeredgewidth(lw) # tick line width

    for line in (ax.xaxis.get_ticklines(minor=True) + ax.yaxis.get_ticklines(minor=True)):
        line.set_markersize(tick_length/2) # tick length
        line.set_markeredgewidth(lw/2) # tick line width

################################# Differential Expression ################################

def _run(x):
    "Generic run."
    res = subprocess.run(x.split(' '), stdout=PIPE, stderr=PIPE, universal_newlines=True)
    res.stdout = res.stdout.strip().split('\n')
    res.stderr = res.stderr.strip().split('\n')
    return res

def write_log(fout,res):
    """Write output logs."""
    with open(fout+'.out', "w") as f:
        for l in res.stdout:
            f.write(l+'\n')

    with open(fout+'.err', "w") as f:
        for l in res.stderr:
            f.write(l+'\n')

def run_differential_expression(formats):
    """
    Run differential expression.
    formats: parameter dictionary
    - input: samples x categories membership matrix
    - covar: samples x covariates matrix
    - run_name
    - clust
    - prot_dir
    - rna_counts
    - out_limma
    - sva: boolean
    """
    if formats['clust'] != 'cohort':
        cmd = "Rscript {}/run_limma.R".format("../scripts")
        cmd += " --covar {}".format(formats['covar'])
    else:
        cmd = "Rscript {}/run_limma_no_covar.R".format("../scripts")
    if formats['sva']:
        cmd += ' --sva TRUE'
    else:
        cmd += ' --sva FALSE'

    cmd += ' --labels {}'.format(formats['input'])\
        +' --label_id {}'.format(formats['clust'])\
        +' --feature_maps {}'.format(formats['prot_maps'])\
        +' --proteome {}/proteome_X.tsv'.format(formats['prot_dir'])\
        +' --phosphoproteome {}/phosphoproteome_X.tsv'.format(formats['prot_dir'])\
        +' --transcriptome {}'.format(formats['rna_counts'])\
        +' --output {}'.format(formats['out_limma'])
    print(cmd)

    res = _run(cmd)
    write_log(os.path.join(formats['out_limma'], "de.log"), res)

    cmd = "python3 {}/postprocess_limma_de.py".format("../scripts")
    cmd += ' -i {}'.format(formats['out_limma'])\
        +' -f {}'.format(formats['prot_maps'])\
        +' -o {}'.format(formats['out_limma'])

    res = _run(cmd)
    write_log(os.path.join(formats['out_limma'], "de.postprocess.log"), res)

def run_differential_expression_cohort(formats):
    """
    Run differential expression.
    """
    cmd = "Rscript {}/run_limma_no_covar.R".format("../scripts")
    if formats['sva']:
        cmd += ' --sva TRUE'
    else:
        cmd += ' --sva FALSE'

    if 'prot_path' not in formats:
        formats['prot_path'] = "proteome_X.tsv"
    if 'phospho_path' not in formats:
        formats['phospho_path'] = "phosphoproteome_X.tsv"

    cmd += ' --labels {}'.format(formats['input'])\
        +' --label_id {}'.format(formats['clust'])\
        +' --feature_maps {}'.format(formats['prot_maps'])\
        +' --proteome {}/{}'.format(formats['prot_dir'],formats['prot_path'])\
        +' --phosphoproteome {}/{}'.format(formats['prot_dir'],formats['phospho_path'])\
        +' --transcriptome {}'.format(formats['rna_counts'])\
        +' --output {}'.format(formats['out_limma'])
    if 'acetyl_path' in formats:
        cmd += ' --acetylome {}/{}'.format(formats['prot_dir'],formats['acetyl_path'])
    print(cmd)

    res = _run(cmd)
    write_log(os.path.join(formats['out_limma'], "de.log"), res)

    cmd = "python3 {}/postprocess_limma_de.py".format("../scripts")
    cmd += ' -i {}'.format(formats['out_limma'])\
        +' -f {}'.format(formats['prot_maps'])\
        +' -o {}'.format(formats['out_limma'])

    res = _run(cmd)
    write_log(os.path.join(formats['out_limma'], "de.postprocess.log"), res)

################################### GSEA #####################################
def write_gct(X, obs, var, outfile):
    """
    Write GCT File.
    """
    hd = list((X.shape[0], X.shape[1], var.shape[1], obs.shape[1]))
    hd = [str(x) for x in hd]
    hd = '\t'.join(hd) + '\n'

    # Build GCT File
    gct_df = pd.concat([var, X],1)
    gct_df = pd.concat((obs.T,gct_df)).loc[:,gct_df.columns]

    # Create file
    with open(outfile, 'w') as f:
        f.write("#1.3\n")
        f.write(hd)
        gct_df.to_csv(f, sep='\t')

def gen_phospho_ptmgsea(run_name, map_df, outfile, weights='gsea_rank'):
    """
    Generate input file for run name.
    """
    from ast import literal_eval

    de_df = pd.read_csv(os.path.join('/home/yakiyama/DE_results', run_name, 'full_diffexp_results.tsv'), sep='\t', index_col=0)
    de_df = de_df[de_df['feature']=='phosphoproteome']

    #de_df = de_df[de_df['qval']<0.1]

    de_df = de_df.loc[:,['gene_name', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B','qval', 'id', 'gsea_rank']]
    de_df = de_df.join(map_df)

    # Phospho
    phosph_var = de_df[['geneSymbol','id.description','accession_number','protein_mw','variableSites','sequence','sequenceVML','VMsiteFlanks']].drop_duplicates()
    phosph_var['VMsiteFlanks'] = phosph_var['VMsiteFlanks'].apply(literal_eval)
    phosph_var = phosph_var.explode('VMsiteFlanks')
    phosph_var['ptmGSEA'] = phosph_var['VMsiteFlanks'].str.upper()+'-p'
    phosph_var = phosph_var.reset_index().drop_duplicates(subset=['index','VMsiteFlanks']).set_index('index')

    phosph_X = de_df.reset_index()[['index',weights,'id']].pivot(index='index',columns='id')
    phosph_X.columns = phosph_X.columns.droplevel(0)
    phosph_X = phosph_X.join(phosph_var['ptmGSEA']).set_index("ptmGSEA")

    phosph_var = phosph_var.reset_index().set_index("ptmGSEA")
    phosph_obs = pd.DataFrame(phosph_X.columns, columns=['S']).set_index("S")
    phosph_obs['run_name'] = run_name
    print(phosph_X)

    # Write GCT
    write_gct(phosph_X, phosph_obs, phosph_var, outfile)

def runPTMSEA(ssGSEA_PATH, gct_input, output_dir):
    cmd = ''
    cmd = 'Rscript ' + ssGSEA_PATH + "ssgsea-cli.R"
    assert(gct_input)
    db_path = ssGSEA_PATH + "db/ptmsigdb/ptm.sig.db.all.flanking.human.v1.9.0.gmt"
    out_prefix = output_dir + '/out'
    cmd += f" -i {gct_input} -o {out_prefix} --db {db_path}"
    print(cmd)
    res = _run(cmd)
    write_log(os.path.join(output_dir, "ptm_sea.log"), res)

def t_test(mat, group_s, equal_var=False):
    """
    t-test compatible w/ Missing Values
    ---------------------
    Args:
        * mat: pd.DataFrame (genes x samples)
        * group_s: series of groupings
    """
    from scipy import stats
    from statsmodels.stats.multitest import multipletests
    from tqdm import tqdm
    mat = mat[group_s.index]
    def _collapser(x, index, columns, name):
        _df = pd.DataFrame(x, index=index, columns=columns).reset_index()
        _id = _df.columns[0]
        return pd.melt(
            pd.DataFrame(x, index=index, columns=columns).reset_index(),
            id_vars=_id,
            ).set_index(_id).rename(columns={'variable':group_s.name,'value':name})
    groups = np.array(group_s)
    X = mat.values
    n_groups = np.unique(groups).shape[0]
    n_genes = X.shape[0]
    # Init np.arrays
    t_stat = np.zeros((n_genes, n_groups))
    pval = np.zeros((n_genes, n_groups))
    pval_adj = np.zeros((n_genes, n_groups))
    x_in = np.zeros((n_genes, n_groups))
    x_out = np.zeros((n_genes, n_groups))
    for idx,group in tqdm(enumerate(np.unique(groups)), total=n_groups):
        mask = groups==group
        if sum(mask) > 1:
            X_in = X[:,mask]
            X_out = X[:,~mask]
            t_stat[:,idx], pval[:,idx] = stats.ttest_ind(X_in, X_out, axis=1, equal_var=equal_var, nan_policy='omit')
            avail_pval = ~np.isnan(pval[:,idx])
            _,pval_adj[avail_pval,idx],_,_ = multipletests(
                pval[avail_pval,idx],
                alpha=0.05,
                method='fdr_bh',
                is_sorted=False,
                returnsorted=False
            )
            x_in[:,idx] = np.nanmean(X_in,1)
            x_out[:,idx] = np.nanmean(X_out,1)
    # Collapse to dataframe
    de_df = pd.concat([
                _collapser(x_in, mat.index, np.unique(groups), 'x_in'),
                _collapser(x_out, mat.index, np.unique(groups), 'x_out')['x_out'],
                _collapser(t_stat, mat.index, np.unique(groups), 't')['t'],
                _collapser(pval, mat.index, np.unique(groups), 'pval')['pval'],
                _collapser(pval_adj, mat.index, np.unique(groups), 'pval_adj')['pval_adj'],
            ],1)
    # Fold-change
    de_df['diff'] = de_df['x_in'] - de_df['x_out']
    # Signed FC * -log10(qval)
    de_df['gsea_rank'] = de_df['diff'] * -np.log10(de_df['pval_adj'])
    # Fraction Missing
    de_df = de_df.join(pd.DataFrame(mat.isna().sum(1) / mat.shape[1], columns=['frac_missing']))
    return de_df

def prepGSEA(de_df, feature, group, rank_col='gsea_rank'):
    df = de_df[(de_df['feature']==feature) & (de_df['id']==group)].sort_values(
        by=rank_col,ascending=False).set_index('gene_name')[[rank_col]]
    if feature in ['phosphoproteome','phosphoproteome_res','acetylome','acetylome_res']:
        df['abs_rank'] = df[rank_col].map(abs)
        df = df.sort_values('abs_rank',ascending=False).reset_index().drop_duplicates(subset='gene_name',keep='first').set_index('gene_name')[[rank_col]]
        df = df.sort_values(rank_col,ascending=False)
    return(df)

def ptm_pval_fdr(res_df, method='fdr_bh'):
    from statsmodels.stats.multitest import multipletests
    res_df['uniqueIndex'] = res_df.apply(lambda x: x.name + '_res' if '_res' in x['feature'] else x.name, 1)
    res_df['collapsed.adj.P.Val'] = np.nan

    # Apply correction for each feature in the results dataframe
    for feature, feature_df in res_df.groupby('feature'):
        temp_df = pd.DataFrame(columns=['P.Value', 'adj.P.Value', 'allIndices'])
        # For each gene name, look at all corresponding rows (transcript, protein, PTM site), and select the one with the smallest
        # p-value, and multiply by the total number of elements inspected.
        for gene, gene_df in feature_df.groupby('gene_name'):
            mostSig = gene_df['P.Value'].idxmin()
            allIndices_s = gene_df['uniqueIndex']
            # Add maximally significant element to temporary dataframe
            temp_df.loc[mostSig, 'P.Value'] = min(gene_df.loc[mostSig, 'P.Value'] * gene_df.shape[0], 1)
            temp_df.loc[mostSig, 'allIndices'] = allIndices_s.to_list()
        # Run FDR on collapsed p-values
        temp_df['adj.P.Val'] = multipletests(temp_df['P.Value'], method=method, is_sorted=False,
                                             returnsorted=False)[1]
        # Assign collapsed adj.P.Val to original dataframe for each element considered
        temp_dict = temp_df.explode('allIndices').set_index('allIndices').to_dict()['adj.P.Val']
        res_df['collapsed.adj.P.Val'] = res_df.apply(lambda x: temp_dict[x['uniqueIndex']] if x['uniqueIndex'] in temp_dict else
                                                     x['collapsed.adj.P.Val'], 1)
    res_df.drop(columns=['uniqueIndex'], inplace=True)
    return res_df

def loadGoGeneSet(filePath: str):
    # opens GO geneset and stores genes in list
    with open(filePath, "r") as oFile:
        start = 0
        while not start:
            line = oFile.readline().strip()
            if line.startswith('>'):
                start = 1
        lines = oFile.readlines()
        gs_l = [line.strip() for line in lines]
    return gs_l
