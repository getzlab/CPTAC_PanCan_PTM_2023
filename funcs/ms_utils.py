# Author: Yo Akiyama (Dec. 2020)

# Import relevant libraries
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import numpy as np
from sys import stdout
import matplotlib.pyplot as plt
from adjustText import adjust_text
import math
import itertools

# Get dimensions for plot
def get_dims(num_plots):
    rows = int(math.sqrt(num_plots))
    cols = math.ceil(math.sqrt(num_plots))
    if rows*cols < num_plots:
        rows+=1
    return rows, cols

# Construct a Signatures x Proteins correlation, p-value, and q-val matrices
def get_all_corr_p(sig_df, expr_df):
    corr_mat = pd.DataFrame(index=expr_df.index, columns=sig_df.columns)
    p_mat = pd.DataFrame(index=expr_df.index, columns=sig_df.columns)
    adjp_mat = pd.DataFrame(index=expr_df.index, columns=sig_df.columns)
    count = 0
    num_feat = expr_df.shape[0]

    # Double loop to fill Correlation and P-val matrices
    for sig in sig_df.columns:
        sig_s = sig_df[sig_df[sig]>0][sig]
        stdout.write('---------Correlation and P-val Matrices: Evaluating {}---------\n'.format(sig))
        for feat in expr_df.index:
            # Subset for shared patients
            feat_s = expr_df.loc[feat, sig_s.index].dropna()
            feat_s = feat_s.loc[np.intersect1d(feat_s.index, sig_s.index)]
            temp_sig_s = sig_s.loc[feat_s.index]
            count += 1
            stdout.write('\r {}/{} proteins'.format(count,num_feat))
            corr_mat.loc[feat,sig], p_mat.loc[feat,sig] = stats.spearmanr(feat_s, temp_sig_s)
        _,adjp,_,_ = multipletests(p_mat[sig], alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)
        # Add missing features and mark with NaN
        # adjp_s = pd.Series(adjp, index=feat_s.index)
        # notIncluded = [x for x in expr_df.index if x not in adjp.index]
        # a
        # djp = pd.concat([adjp, pd.Series(index=notIncluded)])
        adjp_mat[sig] = adjp
        count = 0
        stdout.write('\n')

    return corr_mat, p_mat, adjp_mat

def plot_volcano(corr_vec, q_vec, ax, sig, protein_map=None, numlab=20, isCorr=True, plotGenes=None, plotAll=False, sigOnly=False):
    # Prep plotting dataframe with color labels
    combined_df = pd.DataFrame({
        'corr': corr_vec.astype(float),
        'q-val': q_vec.astype(float)})
    def add_color(row):
        if row['corr'] > 0 and row['q-val'] <= 0.1:
            return 'red'
        elif row['corr'] < 0 and row['q-val'] <= 0.1:
            return 'blue'
        else:
            return 'grey'
    combined_df['color'] = combined_df.apply(add_color,axis=1)

    if protein_map is not None:
        combined_df.index = combined_df.index.map(lambda x: protein_map.loc[x,'geneSymbol'])
    
    # Select set of points to be labeled
    upreg = combined_df[(combined_df.loc[:,'corr'] > 0) & (combined_df.loc[:,'q-val'] <= 0.1)].sort_values(by='q-val')
    downreg = combined_df[(combined_df.loc[:,'corr'] < 0) & (combined_df.loc[:,'q-val'] <= 0.1)].sort_values(by='q-val')
    if upreg.shape[0] > numlab:
        upreg = upreg.iloc[:numlab]
    if downreg.shape[0] > numlab:
        downreg = downreg.iloc[:numlab]

    if ((plotGenes is not None) and (plotAll)):
        subset_df = combined_df.loc[np.intersect1d(plotGenes,combined_df.index)]
        if sigOnly:
            subset_df = subset_df[subset_df['q-val']<=0.1]
        add_upreg = subset_df[subset_df['corr'] > 0]
        add_downreg = subset_df[subset_df['corr'] < 0]
        upreg = pd.concat([upreg, add_upreg]).drop_duplicates(keep='first').sort_values('q-val')
        downreg = pd.concat([downreg, add_downreg]).drop_duplicates(keep='first').sort_values('q-val')
    elif ((plotGenes is not None) and (not plotAll)):
        subset_df = combined_df.loc[np.intersect1d(plotGenes,combined_df.index)]
        if sigOnly:
            subset_df = subset_df[subset_df['q-val']<=0.1]
        upreg =  subset_df[subset_df['corr'] > 0].sort_values('q-val')
        downreg =  subset_df[subset_df['corr'] < 0].sort_values('q-val')
        
    ###### Plotting Block #####
    ax.scatter(combined_df['corr'], -1*np.log(combined_df['q-val']), color=combined_df['color'])
    if isCorr:
        ax.set_xlim(left=-1, right=1)
    #ax[row,col].set_ylim(bottom=-0.01, top=max(-1np.log

    # Add labels for top 20 up/downregulated genes and reduce overlap of labels
    texts = []
    for corr, q, name in zip(upreg['corr'],upreg['q-val'],upreg.index):
        prot_sym = name
        texts.append(ax.text(corr, -1*np.log(q), prot_sym))        
    for corr,q,name in zip(downreg['corr'],downreg['q-val'], downreg.index):
        prot_sym = name
        texts.append(ax.text(corr, -1*np.log(q), prot_sym))
    adjust_text(texts,arrowprops=dict(arrowstyle="->", color='r', lw=0.5),ax=ax)
    # Title with signature 
    ax.set_title(sig)

    return 0

def proteome_mutation_analysis(
        expr_df: pd.DataFrame,
        sig_df: pd.DataFrame,
        feature_type: str,
        protein_map: pd.DataFrame = None,
        out_path: str = None,
        numlab: int = 20,
        fig_height: int = 16,
        fig_width: int = 16,
        axes_fontsize: int = 10,
        title_fontsize: int = 16
):
    # Make sure expression and signature dataframes contain the same samples
    try:
        sig_df.index = expr_df.index
    except:
        raise Exception("Expression and signature dataframes contain different sample names.")
    
    NUMSIGS = sig_df.shape[1]
    stdout.write('{} Active Signatures\n'.format(NUMSIGS))

    # Evaluate correlations and significance
    corr_mat, p_mat, q_mat = get_all_corr_p(sig_df,expr_df)
    
    # Initialize figure with N plots (N = # of Signatures)
    NUMROWS, NUMCOLS = get_dims(NUMSIGS)
    fig, ax = plt.subplots(NUMROWS, NUMCOLS, figsize=(fig_height,fig_width))
    row, col = 0, 0

    # Iterate through signatures to create each volcano plot
    # If 1 signature, no loop
    if NUMSIGS == 1:
        sig = corr_mat.index[0]
        stdout.write('----------------------------Plotting {}----------------------------\n'.format(sig))
        plot_volcano(corr_mat.loc[sig], q_mat.loc[sig], ax, sig, protein_map, numlab)
    # If 2 signatures, 1d axis array
    elif NUMSIGS == 2:
        for sig in corr_mat.index:
            stdout.write('----------------------------Plotting {}----------------------------\n'.format(sig))
            plot_volcano(corr_mat.loc[sig], q_mat.loc[sig], ax[col], sig, protein_map, numlab)
            col+=1
    # For more than 2 signatures, index 2d axis array
    else:
        for sig in corr_mat.index:
            stdout.write('----------------------------Plotting {}----------------------------\n'.format(sig))
            plot_volcano(corr_mat.loc[sig], q_mat.loc[sig], ax[row,col], sig, protein_map, numlab)
            if col//(NUMCOLS-1) == 1:
                row += 1
            col = (col+1) % NUMCOLS

        while col < NUMCOLS-1:
            fig.delaxes(ax[row,col])
            col += 1
        
    # Label figure with titles
    fig.suptitle('Signature-{} Correlations'.format(feature_type), fontsize=title_fontsize,y=0.95)
    fig.text(0.5, 0.1, 'Spearman rho', ha='center', va='center', fontsize=axes_fontsize)
    fig.text(0.1, 0.5, '-log(Adj. P-value)', ha='center', va='center', rotation='vertical',fontsize=axes_fontsize)
    plt.subplots_adjust(wspace=0.1, hspace=0.1)

    # Save files if output path is specified
    if out_path:
        corr_mat.to_csv('{}/corr_table.tsv'.format(out_path),sep='\t')
        p_mat.to_csv('{}/pval_table.tsv'.format(out_path),sep='\t')
        q_mat.to_csv('{}/adj_pval_table.tsv'.format(out_path),sep='\t')
        fig.savefig('{}/volcano_plot.tsv'.format(out_path))
    return corr_mat,p_mat,q_mat,fig

def prep_plotting_df(H, color_scheme, cohort_df, signature):
    cohort_dict = dict(zip(cohort_df['sample'],cohort_df['cohort']))

    H_norm = H.div(H.sum(1), axis=0)
    sig_norm = H_norm[H_norm[signature] > 0].sort_values(by=signature)
    sig_raw = H[H[signature] > 0].sort_values(by=signature)

    sig_plot = pd.concat([sig_raw[signature].rename('raw'), 
                          sig_norm[signature].rename('norm')], axis=1)
    sig_plot['cohort'] = sig_plot.index.map(cohort_dict)
    sig_plot['color'] = sig_plot['cohort'].map(color_scheme)
    return sig_plot

def plot_norm_raw_sigs(sig_plot, raw_thresh=None, norm_thresh=None, circle=None, sharex=True, sharey=True,
                       height=12,width=12,title=None, show_raw_thresh=False, show_norm_thresh=False, plot_all=True,
                       xlabel=None, ylabel=None, xlabel_coords=(0.5, 0.05), ylabel_coords=(0.025, 0.2), nrows=None, ncols=None,
                       fig_suptitle_ycoord = 0.95, fig_suptitle_size=20, ylabel_size=16, xlabel_size=16, show_legend=True, alpha=0.5):
    assert(((nrows is None) and (ncols is None)) or ((nrows is not None) and (ncols is not None)), "nrows and ncols must be both None or not None")
    nCohorts = sig_plot['cohort'].nunique()
    cohorts = sig_plot['cohort'].unique()
    if ((nrows is not None) and (ncols is not None)):
        assert(nrows*ncols > nCohorts, "nrows * ncols must be greater than the number of cohorts")
    if plot_all:
        if ((nrows is None) and (ncols is None)):
            nrows, ncols = get_dims(nCohorts + 1)
        subplot_loc = list(itertools.product(range(nrows), range(ncols)))[1:nCohorts + 1]
    else:
        if ((nrows is None) and	(ncols is None)):
            nrows, ncols = get_dims(nCohorts)
        subplot_loc = list(itertools.product(range(nrows), range(ncols)))[:nCohorts]
    if nrows == 1:
        subplot_loc = [i[1] for i in subplot_loc]
    subplot_dict = {cohorts[i]: subplot_loc[i] for i in range(nCohorts)}

    fig, ax = plt.subplots(nrows, ncols, figsize=(width, height), sharex=sharex, sharey=sharey)
    groups = sig_plot.groupby('cohort')
    labels = []
    subplots = []
    for name, group in groups:
        labels.append(name)
        if circle and plot_all:
            all_co = ax[0, 0].scatter(group['norm'], group['raw'], c=group['color'],
                                      edgecolors=group[circle].map(lambda x: 'black' if x==True else 'none'), alpha=alpha)
            subplots.append(all_co)
        elif plot_all:
            all_co = ax[0, 0].scatter(group['norm'], group['raw'], c=group['color'], alpha=alpha)
            subplots.append(all_co)
        if circle:
            subplot_group = ax[subplot_dict[name]].scatter(group['norm'], group['raw'], c=group['color'], label=name,
                                                           edgecolors=group[circle].map(lambda x: 'black' if x==True else 'none'),
                                                           alpha=alpha)
        else:
            subplot_group = ax[subplot_dict[name]].scatter(group['norm'], group['raw'], c=group['color'], label=name,
                                                           alpha=alpha)
        subplots.append(subplot_group)
        if show_raw_thresh:
            ax[subplot_dict[name]].axhline(y=raw_thresh, linestyle='--', label='_nolegend_')
        if show_norm_thresh:
            ax[subplot_dict[name]].axvline(x=norm_thresh, linestyle='--', label='_nolegend_')

    if show_raw_thresh and plot_all:
        ax[0, 0].axhline(y=raw_thresh, linestyle='--', label='_nolegend_')
    if show_norm_thresh and plot_all:
        ax[0, 0].axvline(x=norm_thresh, linestyle='--', label='_nolegend_')
    if show_legend:
        leg = fig.legend(subplots, labels=labels, loc="upper left", borderaxespad=0.1, title="Cohorts",
                         bbox_to_anchor=(0.15, 0.85), prop={'size':12})
        leg.get_title().set_fontsize(14)
    if xlabel is None:
        xlabel = 'Fraction of mutations attributed to signature'
    if ylabel is None:
        ylabel = 'Number of mutations attributed to signature'
    fig.text(xlabel_coords[0], xlabel_coords[1], xlabel, ha='center', fontsize=xlabel_size)
    fig.text(ylabel_coords[0], ylabel_coords[1], ylabel, ha='center', rotation='vertical', fontsize=ylabel_size)
    fig.suptitle(title if title else 'Signature Attribution Plot', y=fig_suptitle_ycoord, size=fig_suptitle_size)

    return(fig, ax)

def print_quickstats(df, cohort_df, raw_thresh, norm_thresh):
    cohorts = df['cohort'].unique()
    for group in cohorts:
        proportion_raw = df[(df['cohort'] == group) & (df['raw'] >= raw_thresh)].shape[0] / cohort_df[cohort_df['cohort'] == group].shape[0]
        proportion_frac = df[(df['cohort'] == group) & (df['norm'] >= norm_thresh)].shape[0] / cohort_df[cohort_df['cohort'] == group].shape[0]
        print("{}:\n\tPasses raw threshold ({} mutations): {}\n".format(group, raw_thresh, proportion_raw))
        print("\tPasses fractional threshold ({} of mutations): {}\n".format(group, norm_thresh, proportion_frac))

