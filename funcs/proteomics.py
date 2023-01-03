import sys
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from typing import Union

import nmf_utilities as nmu
from utils import get_pcs

# ----------------------------------
# Proessing Proteomics
# ----------------------------------
class GCT(object):
    """
    GCT python class.
    ------------------
    Class for reading and storing GCT data.
    """
    def __init__(self, filepath: str, rm_var_whitespace: bool = True, parse_fn: bool = True):
        """
        Args:
            * filepath: filename
        """
        with open(filepath) as f:
            self.n_var, self.n_obs, self.ny, self.nx = map(int,f.readlines()[1].strip().split("\t"))

        self.filename = os.path.basename(filepath)
        self.X = pd.read_csv(filepath, sep='\t', index_col=0, skiprows=2).iloc[self.nx:,self.ny:]
        self.var = pd.read_csv(filepath, sep='\t', index_col=0, skiprows=2).iloc[self.nx:,:self.ny]
        self.obs = pd.read_csv(filepath, sep='\t', index_col=0, skiprows=2).iloc[:self.nx,self.ny:].T

        if rm_var_whitespace:
            self.X.index = self.X.index.str.replace(" ","")
            self.var.index = self.var.index.str.replace(" ","")

        self.var_names = np.array(self.var.index)
        self.obs_names = np.array(self.obs.index)

        if parse_fn:
            self.cohort, self.assay = os.path.basename(filepath).split('-')[1:3]
            self.cohort = self.cohort.upper()

            self.var.loc[:,'feature'] = self.assay
            self.obs.loc[:,'cohort'] = self.cohort
        else:
            self.cohort, self.assay = None, None

    def __str__(self):
        if self.cohort is None:
            return "GCT | ({} feats x {} obs)".format(self.n_var, self.n_obs)
        else:
            return "GCT | {} - {} ({} feats x {} obs)".format(self.cohort, self.assay, self.n_var, self.n_obs)

    def subset(self, obs=None, var=None):
        """Subset indices."""
        if obs is not None:
            self.X = self.X.loc[:,obs]
            self.obs = self.obs.loc[obs,:]
            self.obs_names = np.array(self.obs.index)
            self.n_obs = self.obs_names.shape[0]

        if var is not None:
            self.X = self.X.loc[var,:]
            self.var = self.var.loc[var,:]
            self.var_names = np.array(self.var.index)
            self.n_var = self.var_names.shape[0]

    def rename(self, obs_d=None, var_d=None):
        """Rename indices."""
        if obs_d is not None:
            self.X = self.X.rename(columns=obs_d)
            self.obs = self.obs.rename(index=obs_d)
            self.obs_names = np.array(self.obs.index)
            self.n_obs = self.obs_names.shape[0]

        if var_d is not None:
            self.X = self.X.rename(index=var_d)
            self.var = self.var.rename(index=var_d)
            self.var_names = np.array(self.var.index)
            self.n_var = self.var_names.shape[0]

def collapse_prot_to_gene(g: GCT):
    """
    Collapse Protein to Gene-Centric.
    """
    import re

    for c in ('numSpectraProteinObserved','subgroupNum'):
        assert c in g.var.columns, "{} missing in GCT file.".format(c)

    # Cast columns
    _var = g.var.copy()
    _var['numSpectraProteinObserved'] = _var['numSpectraProteinObserved'].astype(int)
    _var['proteinGroup'] = _var['subgroupNum'].apply(lambda x: str(x).replace(".","_").split("_")[0]).astype(int)
    _var['proteinRank'] = _var['subgroupNum'].apply(lambda x: str(x).replace(".","_").split("_")[1]).astype(int)
    _var = _var.sort_values(['proteinGroup','proteinRank'])

    # Select GeneSymbol
    _p = _var.shape[0]
    _var = _var.drop_duplicates(subset=['geneSymbol','proteinGroup'])

    # Unique GeneSymbol
    _var = _var.sort_values('numSpectraProteinObserved', ascending=False).drop_duplicates('geneSymbol')

    print("      * {} / {} protein rows filtered".format(_p-_var.shape[0], _p))
    g.var['id'] = g.var.index
    g.subset(var=_var.index)
    g.rename(var_d=_var['geneSymbol'].to_dict())

def normal_name(x):
    if isinstance(x,str):
        if "NAT" in x or "Normal" in x or "normal" in x:
            return "Normal"
        elif "Tumor" in x:
            return "Tumor"
    else:
        return None

def _melt_ptmgsea(df, col, group_names, fmt=''):
    """
    Return melted.
    """
    import re

    df = df.filter([col+"."+g for g in group_names]).copy().reset_index().melt(id_vars=['id'])
    df['variable'] = df['variable'].apply(lambda x: x.split(col+'.')[1])
    return df.rename(columns={'value':col, 'variable':'signature','id':'pathway'})

def gen_ptmgsea_df(gct, **kwargs):
    """
    Generate full PTM-GSEA table.
    """
    s = ['Signature.set.overlap.percent','Signature.set.overlap','pvalue','fdr.pvalue']

    X = gct.X.copy()
    group_names = X.columns
    X.columns = 'NES.' + X.columns

    ptm_gsea_df = _melt_ptmgsea(
        gct.var,s[0],group_names, **kwargs
    ).merge(
        _melt_ptmgsea(gct.var,s[1],group_names, **kwargs)
    ).merge(
        _melt_ptmgsea(gct.var,s[2],group_names, **kwargs)
    ).merge(
        _melt_ptmgsea(gct.var,s[3],group_names, **kwargs)
    ).merge(
        _melt_ptmgsea(X,"NES",group_names, **kwargs)
    )

    return ptm_gsea_df.rename(columns={'signature':'id'})

def parse_participant_and_type(gct,assay):
    """
    Parse participants and sample type.
    """
    # Add Participant_ID
    gct['MEDUL'][assay].obs['Participant_ID'] = gct['MEDUL'][assay].obs['Sample.ID']
    gct['UCEC'][assay].obs['Participant_ID'] = gct['UCEC'][assay].obs['Sample']
    gct['LSCC'][assay].obs['Participant_ID'] = gct['LSCC'][assay].obs['Participant']
    gct['LUAD'][assay].obs['Participant_ID'] = gct['LUAD'][assay].obs['Participant']
    gct['BRCA'][assay].obs['Participant_ID'] = gct['BRCA'][assay].obs['Sample.ID'].apply(lambda x: x.replace('X','').split('.')[0])
    gct['GBM'][assay].obs['Participant_ID'] = gct['GBM'][assay].obs['Sample.ID'].str.replace('.','-')

    if assay is not 'acetylome':
        gct['HNSCC'][assay].obs['Participant_ID'] = gct['HNSCC'][assay].obs['Participant']
        gct['PDAC'][assay].obs['Participant_ID'] = gct['PDAC'][assay].obs['X1211.CPTAC_Participant_ID'].str.replace('.','-')
        gct['CCRCC'][assay].obs['Participant_ID'] = gct['CCRCC'][assay].obs['Case_ID']
        gct['OV'][assay].obs['Participant_ID'] = gct['OV'][assay].obs['CPTAC_Case_ID']
        gct['COAD'][assay].obs['Participant_ID'] = gct['COAD'][assay].obs['Sample.ID.noNC']

    # Fix Type
    gct['MEDUL'][assay].obs['Type'] = 'Tumor'

    # Type ID
    gct['MEDUL'][assay].obs['Type_ID'] = gct['MEDUL'][assay].obs['Type']
    gct['UCEC'][assay].obs['Type_ID'] = gct['UCEC'][assay].obs['Type']
    gct['LSCC'][assay].obs['Type_ID'] = gct['LSCC'][assay].obs['Type'].apply(normal_name)
    gct['LUAD'][assay].obs['Type_ID'] = gct['LUAD'][assay].obs['Type'].apply(normal_name)
    gct['BRCA'][assay].obs['Type_ID'] = gct['BRCA'][assay].obs['Type'].apply(normal_name)
    gct['GBM'][assay].obs['Type_ID'] = gct['GBM'][assay].obs['Type'].apply(normal_name)

    if assay is not 'acetylome':
        gct['HNSCC'][assay].obs['Type_ID'] = gct['HNSCC'][assay].obs['Type'].apply(normal_name)
        gct['PDAC'][assay].obs['Type_ID'] = gct['PDAC'][assay].obs['Type'].apply(normal_name)
        gct['CCRCC'][assay].obs['Type_ID'] = gct['CCRCC'][assay].obs['Type']
        gct['OV'][assay].obs['Type_ID'] = gct['OV'][assay].obs['Type']
        gct['COAD'][assay].obs['Type_ID'] = gct['COAD'][assay].obs['Type'].apply(normal_name)

def compile_feature_counts(gct, assay):
    """
    Compile feature counts.
    -----------------------
    Args:
        * gct: dictionary of gct files
    """
    _df_s = list()
    for k in gct.keys():
        if assay in gct[k].keys():
            _df = pd.DataFrame(gct[k][assay].var['geneSymbol'])
            _df['{}_frac'.format(k)] = 1-gct[k][assay].X.isna().sum(1) / gct[k][assay].n_obs
            _df['{}_n_samples'.format(k)] = gct[k][assay].n_obs
            _df_s.append(_df)

    res_df = pd.concat(_df_s, axis=1).groupby(level=0, axis=1).first()

    cols = res_df.columns.to_list()
    cols.insert(0, cols.pop(cols.index('geneSymbol')))
    res_df = res_df.reindex(columns=cols)

    return res_df

def subset_pdf(df, tn, tag):
    """
    Subset protein df.
    """
    return df[df['Type_ID']==tn].reset_index().set_index("Participant_ID")[['index']].rename(columns={'index':'{}_Sample_ID'.format(tag)})

def get_pca_umap(X, n_components=25, normalize=False, **kwargs):
    """
    Get PCA UMAP.
    """
    import umap

    P_df, pca, feats = get_pcs(X, normalize=normalize, return_genes=True, n_components=n_components)

    fit = umap.UMAP(**kwargs)
    H_umap = fit.fit_transform(P_df)
    P_df['umap1'] = H_umap[:,0]
    P_df['umap2'] = H_umap[:,1]

    return P_df, pca, feats

def process_ptm_var(pmap_df):
    """
    Process PTM var.
    """
    pmap_df['protein_mw'] = pmap_df['protein_mw'].astype(float)
    pmap_df['variableSites'] = pmap_df['variableSites'].str.strip()
    pmap_df['sequence'] = pmap_df['sequence'].str.strip()
    pmap_df['VMsiteFlanks'] = pmap_df['VMsiteFlanks'].str.strip()
    pmap_df['sequenceVML'] = pmap_df['sequenceVML'].str.strip()

    pmap_df = pmap_df.reset_index().drop_duplicates().set_index('id')

    pmap_seq_df = pmap_df[['variableSites']].groupby(level=0).agg({
        'variableSites': lambda x: list(x)
    }).join(
        pmap_df[['sequence']].groupby(level=0).agg({
            'sequence': lambda x: list(x)
        })
    ).join(
        pmap_df[['VMsiteFlanks']].dropna().groupby(level=0).agg({
            'VMsiteFlanks': lambda x: list(x)
        })
    ).join(
        pmap_df[['sequenceVML']].dropna().groupby(level=0).agg({
            'sequenceVML': lambda x: list(x)
        })
    )

    pmap_df = pmap_df.reset_index().loc[:,[
        'id',
        'id.description',
        'geneSymbol',
        'protein_mw',
        'accession_number',
        'feature'
    ]].drop_duplicates().set_index('id')

    return pmap_df.join(pmap_seq_df)

def assign_ptm_type(x: str):
        if x.startswith("ENSG"):
                return "RNA"
        elif len(x.split('.')[-1]) < 4:
                return "Protein"
        else:
                if "K" in x:
                        return "Acetylation"
                else:
                        return "Phosphorylation"

def dict_matrix_loader(path: str):
        """
        Load matrices into dictionary.
        -----------------------
        Load all PTM / protein features.
        """
        X = {}

        X['acetylome'] = pd.read_parquet(os.path.join(path, "acetylome_X.parquet")).T
        X['proteome'] = pd.read_parquet(os.path.join(path, "proteome_X.parquet")).T
        X['phosphoproteome'] = pd.read_parquet(os.path.join(path, "phosphoproteome_X.parquet")).T

        return X

def map_to_genes(X_i: pd.DataFrame, map_s: pd.Series):
        """
        Map matrix to gene names.
        ------------------------
        Args:
                * X_i: pd.DataFrame (n samples x protein sites)
                * map_s: pd.Series of mapping (i.e. gene_name)

        Returns:
                * pd.DataFrame with mapped feature names
        """
        X = X_i.copy().T
        X['idx'] = np.arange(X.shape[0])
        return X.join(map_s).set_index(map_s.name).sort_values("idx").drop(columns=['idx']).T

def get_enriched_genes(signatures: pd.DataFrame, map_s: pd.Series, diff_thresh: float = 0.5, max_norm: float = 0.5):
        """
        Get Enriched Genes
        ---------------
        Maps proteomics sites to gene names.
        Returns these for upregulated and downregulated genes.
        """
        signatures_p = signatures.loc[[x for x in signatures.index if not x.endswith("_n")]].iloc[:,-5:]
        signatures_p = signatures_p[(signatures_p['diff'] > diff_thresh) & (signatures_p['max_norm'] > max_norm)]

        signatures_n = signatures.loc[[x for x in signatures.index if x.endswith("_n")]].iloc[:,-5:]
        signatures_n = signatures_n[(signatures_n['diff'] > diff_thresh) & (signatures_n['max_norm'] > max_norm)]

        signatures_n.index = signatures_n.reset_index()['index'].apply(lambda x: x.split("_n")[0])

        genes_p = map_to_genes(signatures_p.T, map_s).T
        genes_p = genes_p.loc[genes_p.index.dropna()]
        genes_p['max_id'] = genes_p['max_id'].astype(int)

        genes_n = map_to_genes(signatures_n.T, map_s).T
        genes_n = genes_n.loc[genes_n.index.dropna()]
        genes_n['max_id'] = genes_n['max_id'].astype(int)

        return genes_p, genes_n

def scrape_kw(map_df, tag, map_idx='geneSymbol'):
        """
        Filter mapping dataframe by genes that start
        with tag.
        """
        genes_to_use = [x for x in map_df[map_idx].values if str(x).startswith(tag)]
        return map_df[map_df[map_idx].isin(genes_to_use)]

def scatterplot_by_cohort(X, feat, meta, map_df, map_idx='geneSymbol', tag="", std_thresh=1.5):
        """
        Scatterplot by cohort.
        """
        colors={True:'red', False:'black'}

        x = map_to_genes(X[feat], map_df[map_idx])
        cohort_s = meta[feat]['cohort']
        fig,axes = plt.subplots(1, np.unique(cohort_s).shape[0], figsize=(np.unique(cohort_s).shape[0]*4,4), sharex=True, sharey=True)

        genes_to_keep = []
        for idx,m in enumerate(np.unique(cohort_s)):
                x_m = x.loc[cohort_s[cohort_s==m].index]
                m_df = pd.DataFrame(x_m.mean(0)).rename(columns={0:'mean'})
                m_df['std'] = x_m.std(0)
                m_df['filt'] = m_df.reset_index().apply(lambda x: True if str(x[map_idx]).startswith(tag) and x['std']<std_thresh else False, 1).values
                m_df = m_df.sort_values('filt')

                axes[idx].scatter(m_df['mean'], m_df['std'], alpha=0.5, color=m_df['filt'].apply(lambda x: colors[x]), rasterized=True)
                axes[idx].set_xlabel("$\mu$", fontsize=14)

                if idx==0:
                        axes[idx].set_ylabel("Std. Dev", fontsize=14)

                axes[idx].set_title("{}".format(m), fontsize=16)

                if tag is not "":
                        axes[idx].legend(["norm genes"])

                axes[idx].get_legend().legendHandles[0].set_color('red')
                genes_to_keep.append(np.array(m_df[m_df['filt']==True].index))

        plt.tight_layout()

        if tag is not "":
                genes_to_keep = np.array(list(set.intersection(*[set(x) for x in genes_to_keep])))
                sites_to_keep = X[feat].columns[x.columns.isin(genes_to_keep)]
                plt.suptitle("{} | {} --> {} genes & {} sites.".format(feat, ",".join(tag), genes_to_keep.shape[0], sites_to_keep.shape[0]), fontsize=20, y=1.1)
                return fig, genes_to_keep, sites_to_keep
        else:
                plt.suptitle("{}".format(feat, fontsize=20, y=1.1))
                return fig

        """
        PTM - Protein Comparison
        -----------------
        Args:
                * ptm: pd.DataFrame of PTM
                * prot: pd.DataFrame of Protein expression
        """
        if aggr == 'mean':
                ptm_df = pd.DataFrame(ptm.mean(1)).rename(columns={0:'PTM'})
                prot_df = pd.DataFrame(prot.mean(1)).rename(columns={0:'Proteome'}).reset_index().rename(columns={'id':'protein'})
        elif aggr is None:
                ptm_df = pd.melt(ptm.reset_index(), var_name='sample', id_vars=['id'], value_name='PTM').set_index("id")
                prot_df = pd.melt(prot.reset_index(), id_vars=['id'], var_name='sample', value_name='Proteome').rename(columns={'id':'protein'})

        ptm_df['protein'] = ptm_df.reset_index()['id'].apply(lambda x: "_".join(x.split("_")[:2])).values
        ptm_df = ptm_df.reset_index().merge(prot_df, how='left').set_index("id")

        return ptm_df

def ptm_comparison(ptm, prot, aggr=None, collapsed=False, ptm2Prot_dict=None):
        """
        PTM - Protein Comparison
        -----------------
        Args:
                * ptm: pd.DataFrame of PTM
                * prot: pd.DataFrame of Protein expression
                * aggr: method of aggregation
                        * None - use all datapoints
                        * mean: use mean value (smooths signal, not reccommended)

        Returns:
                * pd.DataFrame
        """
        if aggr == 'mean':
                ptm_df = pd.DataFrame(ptm.mean(1)).rename(columns={0:'PTM'})
                prot_df = pd.DataFrame(prot.mean(1)).rename(columns={0:'Proteome'}).reset_index().rename(columns={'id':'protein'})
        elif aggr is None:
                ptm_df = pd.melt(ptm.reset_index(), var_name='sample', id_vars=['id'], value_name='PTM').set_index("id")
                prot_df = pd.melt(prot.reset_index(), id_vars=['id'], var_name='sample', value_name='Proteome').rename(columns={'id':'protein'})

        if collapsed:
                ptm_df['protein'] = ptm_df.reset_index()['id'].map(ptm2Prot_dict).values
        else:
                ptm_df['protein'] = ptm_df.reset_index()['id'].apply(lambda x: "_".join(x.split("_")[:2])).values
        ptm_df = ptm_df.reset_index().merge(prot_df, how='left').set_index("id")

        return ptm_df

def fit_ptm_prot_ols(ptm, prot, meta=None, aggr=None, ptm2Prot_dict=None, collapsed=False):
        """
        Fit PTM Protein OLS
        -----------------------
        Args:
                * ptm: ptm data
                * prot: proteomics data
                * meta: metadata (cohort)
                * aggr: how to aggregate (default None)

        Return:
                * pd.DataFrame with residuals
                * dict: mapping cohort -> statsmodels fit
        """
        import statsmodels.api as sm
        from tqdm import tqdm

        ptm.index.name = 'id'
        prot.index.name = 'id'

        if meta is None:
                df = ptm_comparison(ptm, prot, aggr=aggr, collapsed=collapsed, ptm2Prot_dict=ptm2Prot_dict)
                df = df.dropna(subset=('PTM','Proteome'))

                # Linear Regression
                mod = sm.OLS(df['PTM'].values[:,np.newaxis], df['Proteome'].values[:,np.newaxis])
                res = mod.fit()
                df['residual'] = res.resid

                print("   * {} / {} sites with matching protein".format(np.unique(df.index).shape[0], ptm.index.shape[0]))
                return df, res
        else:
                regression_dict = {}
                df = list()

                for cohort in tqdm(np.unique(meta), desc='Fitting OLS'):
                        _df = ptm_comparison(ptm[meta[meta==cohort].index], prot[meta[meta==cohort].index], aggr=aggr,
                                             collapsed=collapsed, ptm2Prot_dict=ptm2Prot_dict)
                        _df['cohort'] = cohort
                        _df = _df.dropna(subset=('PTM','Proteome'))

                        # Linear Regression
                        mod = sm.OLS(_df['PTM'].values[:,np.newaxis], _df['Proteome'].values[:,np.newaxis])
                        res = mod.fit()
                        _df['residual'] = res.resid

                        regression_dict[cohort] = res
                        df.append(_df)

                df = pd.concat(df)

                print("   * {} / {} sites with matching protein".format(np.unique(df.index).shape[0], ptm.index.shape[0]))
                return df, regression_dict

def process_Wmat(nmf_file, pmap_df):
    """
    Process W-matrix for pathway enrichment analysis.
    """
    _,_,W,_,_ = nmu.nmf_loader(nmf_file)
    W['n'] = [-1 if x else 1 for x in W.index.str.endswith('_n')]
    W['max'] = W['max']*W['n']
    W['max_norm'] = W['max_norm']*W['n']
    W.index = W.index.str.replace('_n','')
    W = W.join(pmap_df[['geneSymbol','feature']])
    W['qval'] = 0
    W['abs'] = np.abs(W['max_norm'])

    return W

def gen_phospho_wmat_ptmgsea(wmat_file, map_df, outfile, weights='max_norm'):
    """
    Generate input file for run name specifically for W-matrix.
    """
    from ast import literal_eval

    W_df = pd.read_csv(wmat_file, sep='\t', index_col=0)
    W_df = W_df[W_df['feature']=='phosphoproteome'].drop(columns=['feature'])
    W_df = W_df.join(map_df[['id.description','accession_number','protein_mw','variableSites','sequence','sequenceVML','VMsiteFlanks']])

    # Select first
    W_df = W_df.reset_index().sort_values('abs', ascending=False).drop_duplicates(subset='index', keep='first')

    phosph_var = W_df[['index','geneSymbol','id.description','accession_number','protein_mw','variableSites','sequence','sequenceVML','VMsiteFlanks']].drop_duplicates()
    phosph_var['VMsiteFlanks'] = phosph_var['VMsiteFlanks'].apply(literal_eval)
    phosph_var = phosph_var.explode('VMsiteFlanks')
    phosph_var['ptmGSEA'] = phosph_var['VMsiteFlanks'].str.upper()+'-p'
    phosph_var = phosph_var.drop_duplicates(subset=['index']).set_index("index")

    phosph_X = W_df.set_index('index').drop(columns=['max', 'max_id', 'max_norm', 'n', 'geneSymbol',
       'qval', 'abs', 'id.description', 'accession_number', 'protein_mw',
       'variableSites', 'sequence', 'sequenceVML', 'VMsiteFlanks'])
    phosph_X = phosph_X.join(phosph_var['ptmGSEA']).set_index("ptmGSEA")

    phosph_var = phosph_var.reset_index().set_index("ptmGSEA")
    phosph_obs = pd.DataFrame(phosph_X.columns, columns=['S']).set_index("S")
    phosph_obs['run_name'] = wmat_file.split("/")[-2]

    # Write GCT
    write_gct(phosph_X, phosph_obs, phosph_var, outfile)

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

# ----------------------------------
# Sample-Sample Clustering
# ----------------------------------
def ss_linkage_cut(cmatrix, metric='euclidean', method='ward'):
    """
    Get linkage cut from a consensus matrix.
    """
    from scipy.cluster.hierarchy import cut_tree
    import scipy.cluster.hierarchy as shc

    d_linkage = shc.linkage(cmatrix, metric=metric, method=method)
    C = cut_tree(d_linkage)

    return pd.DataFrame(C[:,::-1], index=cmatrix.index)

def ss_linkage_groups(C, lim=40, verbose=False):
    """
    Get groupings for each tree cut.
    """
    res = {}

    def R(X, i, lim=lim):
        """Recurse."""
        if np.unique(X[i]).shape[0]==2:
            if verbose:print("level {} | tot = {}".format(i, X.shape),X.groupby(i).size())
            res[i] = X[[i]]

        if X.shape[0]==1 or X.shape[0]<lim:
            pass
        else:
            for idx in np.unique(X[i]):
                R(X[X[i]==idx],i+1)

    # Run
    R(C,1)

    rename = dict()
    print("  * {} group comparisons".format(len(res.keys())))

    res_df = list()
    for x in res.keys():
        u = np.unique(res[x].iloc[:,0].dropna().copy())
        rename[u[0]] = True
        rename[u[1]] = False
        res_df.append(res[x].applymap(lambda x: rename[x]))
    res_df = pd.concat(res_df,1)
    #return pd.concat([res[x] for x in res.keys()],1)
    return res_df
# ----------------------------------
# Differential expression
# ----------------------------------
import rpy2
from rpy2.robjects.packages import importr
from collections import Iterable

"""
Storey Q-Values - https://github.com/StoreyLab/qvalue
--------------------
Python Wrapper
Author: Francois Aguet
https://github.com/broadinstitute/tensorqtl/blob/master/tensorqtl/rfunc.py
"""
def qvalue(p, lambda_qvalue=None):
        """Wrapper for qvalue::qvalue"""
        qvalue = importr("qvalue")
        rp = rpy2.robjects.vectors.FloatVector(p)
        if lambda_qvalue is None:
                q = qvalue.qvalue(rp)
        else:
                if not isinstance(lambda_qvalue, Iterable):
                        lambda_qvalue = [lambda_qvalue]
                rlambda = rpy2.robjects.vectors.FloatVector(lambda_qvalue)
                q = qvalue.qvalue(rp, **{'lambda':rlambda})
        qval = np.array(q.rx2('qvalues'))
        pi0 = np.array(q.rx2('pi0'))[0]
        return qval, pi0

def t_test(mat, group_s, equal_var=False):
        """
        t-test
        ---------------------
        Args:
                * mat: pd.DataFrame (genes x samples)
                * group_s: series of groupings
        """
        from scipy import stats
        from statsmodels.stats.multitest import multipletests

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
        # qval = np.zeros((n_genes, n_groups))
        x_in = np.zeros((n_genes, n_groups))
        x_out = np.zeros((n_genes, n_groups))

        for idx,group in enumerate(np.unique(groups)):
                mask = groups==group
                if sum(mask) > 1:
                        X_in = X[:,mask]
                        X_out = X[:,~mask]

                        t_stat[:,idx], pval[:,idx] = stats.ttest_ind(X_in, X_out, axis=1, equal_var=equal_var)

                        _,pval_adj[:,idx],_,_ = multipletests(
                                pval[:,idx],
                                alpha=0.05,
                                method='fdr_bh',
                                is_sorted=False,
                                returnsorted=False
                        )

                        #qval[:,idx],_ = qvalue(pval[:,idx])
                        x_in[:,idx] = np.mean(X_in,1)
                        x_out[:,idx] = np.mean(X_out,1)

        # Collapse to dataframe
        de_df = pd.concat([
                                _collapser(x_in, mat.index, np.unique(groups), 'x_in'),
                                _collapser(x_out, mat.index, np.unique(groups), 'x_out')['x_out'],
                                _collapser(t_stat, mat.index, np.unique(groups), 't')['t'],
                                _collapser(pval, mat.index, np.unique(groups), 'pval')['pval'],
                                _collapser(pval_adj, mat.index, np.unique(groups), 'pval_adj')['pval_adj'],
                                #_collapser(qval, mat.index, np.unique(groups), 'qval')['qval']
                        ],1)

        # Fold-change
        de_df['diff'] = de_df['x_in'] - de_df['x_out']

        # Signed FC * -log10(qval)
        de_df['gsea_rank'] = de_df['diff'] * -np.log10(de_df['pval_adj'])

        return de_df

def mannwhitneyu(mat, group_s):
        """
        mannwhitneyu
        ---------------------
        Args:
                * mat: pd.DataFrame (genes x samples)
                * group_s: series of groupings
        """
        from tqdm import tqdm
        from scipy import stats
        from statsmodels.stats.multitest import multipletests
        from sys import stdout

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
        u_stat = np.zeros((n_genes, n_groups))
        pval = np.zeros((n_genes, n_groups))
        pval_adj = np.zeros((n_genes, n_groups))
        qval = np.zeros((n_genes, n_groups))
        x_in = np.zeros((n_genes, n_groups))
        x_out = np.zeros((n_genes, n_groups))

        for idx,group in enumerate(np.unique(groups)):
                stdout.write("\r{} of {}".format(idx+1, n_groups))
                mask = groups==group
                if sum(mask) > 1:
                        X_in = X[:,mask]
                        X_out = X[:,~mask]

                        for gn in range(X_in.shape[0]):
                                u_stat[gn,idx], pval[gn,idx] = stats.mannwhitneyu(X_in[gn], X_out[gn])

                        _,pval_adj[:,idx],_,_ = multipletests(
                                pval[:,idx],
                                alpha=0.05,
                                method='fdr_bh',
                                is_sorted=False,
                                returnsorted=False
                        )

                        try:
                                qval[:,idx],_ = qvalue(fgsea_df['pval'].values)
                        except:
                                try:
                                        qval[:,idx],_ = qvalue(fgsea_df['pval'].values, lambda_qvalue=0.5)
                                except:
                                        qval[:,idx] = None

                        x_in[:,idx] = np.mean(X_in,1)
                        x_out[:,idx] = np.mean(X_out,1)

        # Collapse to dataframe
        de_df = pd.concat([
                                _collapser(x_in, mat.index, np.unique(groups), 'x_in'),
                                _collapser(x_out, mat.index, np.unique(groups), 'x_out')['x_out'],
                                _collapser(u_stat, mat.index, np.unique(groups), 'u')['u'],
                                _collapser(pval, mat.index, np.unique(groups), 'pval')['pval'],
                                _collapser(pval_adj, mat.index, np.unique(groups), 'pval_adj')['pval_adj'],
                                _collapser(qval, mat.index, np.unique(groups), 'qval')['qval']
                        ],1)

        # Fold-change
        de_df['diff'] = de_df['x_in'] - de_df['x_out']

        # Signed FC * -log10(qval)
        de_df['gsea_rank'] = de_df['diff'] * -np.log10(de_df['pval_adj'])

        return de_df

def get_ptm_de_proportions(
        df: pd.DataFrame,
        map_s: Union[None, pd.Series] = None,
        qval_thresh=0.01,
        lfc_thresh=0.5,
        lfc_idx='diff',
        cluster_idx=None,
        norm=True,
        qval_idx='qval'
        ):
        """
        Get PTM DE Proportions
        --------------------------
        Args:
                * de_df: differential expression results
                * map_s: pd.Series of type

        Returns:
                * pd.DataFrame, pd.DataFrame (up-regulated sites, down-regulated sites)
        """
        if cluster_idx is None:
                cluster_idx = df.columns[0]

        if map_s is not None:
                df = df.join(map_s)

        pivot_df_up = nmu.pivot_props(df[(df[qval_idx]<qval_thresh) & (df[lfc_idx]>lfc_thresh)], x='feature', y=cluster_idx, norm=norm)
        #pivot_df_up = pivot_df_up.rename(columns={'Acetylome':'Acetylation', 'Phosphoproteome':'Phosphorylation', 'Proteome':'Protein'})

        pivot_df_down = nmu.pivot_props(df[(df[qval_idx]<qval_thresh) & (df[lfc_idx]<-lfc_thresh)], x='feature', y=cluster_idx, norm=norm)
        #pivot_df_down = pivot_df_down.rename(columns={'Acetylome':'Acetylation', 'Phosphoproteome':'Phosphorylation', 'Proteome':'Protein'})

        return pivot_df_up, pivot_df_down

# ----------------------------------
# fGSEA
# ----------------------------------
def get_gene_set(
        df: pd.DataFrame,
        groupby: str = None,
        how: str = 'mean',
        gene_idx: str ='geneSymbol',
        rank_idx: str ='gsea_rank',
        qval_thresh: float = 0.1,
        qval_idx: str = 'qval'
        ) -> pd.DataFrame:
        """
        Get Gene Set.
        --------------------------
        Args:
                * df: pd.DataFrame of differential expression results
                * groupby: str columns index of groupby
                * how: aggregation method across protein sites
                * gene_idx: gene_symbol column index in df
                * rank_idx: ranking column index in df
                * qval_thresh: q_value threshold for genes to include in the ranked list

        Returns:
                * dictionary mapping each groupby id to ranking pd.DataFrame
        """
        if groupby is None:
                groupby = df.columns[0]

        ranks = dict()

        for group in np.unique(df[groupby]):
                group_df = df[df[groupby]==group]
                group_df = group_df[group_df[qval_idx]<qval_thresh]

                if how=='mean':
                        ranks[group] = group_df[[gene_idx,rank_idx]].groupby(gene_idx).mean()[[rank_idx]].sort_values(rank_idx)
                elif how=='median':
                        ranks[group] = group_df[[gene_idx,rank_idx]].groupby(gene_idx).median()[[rank_idx]].sort_values(rank_idx)
                elif how=='magnitude':
                        def _assign_max(row):
                                if row['m']=='max_abs':
                                        return row['max']
                                elif row['m']=='min_abs':
                                        return row['min']

                        rank_df = group_df[[gene_idx,rank_idx]].groupby(gene_idx).max()[[rank_idx]].rename(columns={rank_idx:'max'}).join(
                                group_df[[gene_idx,rank_idx]].groupby(gene_idx).min()[[rank_idx]].rename(columns={rank_idx:'min'})
                        )

                        rank_df = rank_df.join(rank_df.abs().rename(columns={'max':'max_abs','min':'min_abs'}))
                        rank_df['m'] = rank_df[['max_abs','min_abs']].idxmax(1)
                        rank_df[rank_idx] = rank_df.apply(_assign_max, axis=1)
                        ranks[group] = rank_df[[rank_idx]].sort_values(rank_idx)

        return ranks

def compute_all_enrich(ranks: dict, gmts: list, **kwargs) -> dict:
        """
        Compute All Enrichments
        --------------------------
        Args:
                * ranks: dictionary mapping ID --> pd.DataFrame with index=Genes & column of Ranks
                * gmts: list of GMT files for enrichment
                ** kwargs: passed to rfgsea

        Returns:
                * pd.DataFrame of enrichment results
        """
        from tqdm import tqdm
        enrichment_df = list()

        for k in tqdm(ranks.keys()):
                if sum(ranks[k].isna().values)[0] != ranks[k].shape[0]:
                        e_df = pd.concat([rfgsea(ranks[k], gmt, **kwargs) for gmt in gmts])
                        e_df['id'] = k
                        enrichment_df.append(e_df)

        return pd.concat(enrichment_df)

def rfgsea(r: pd.DataFrame, gmt: str, rank_idx: Union[None, str] = None, nperm=10000, nproc=0, min_size=1) -> pd.DataFrame:
        """
        Wrapper for fgsea::fgsea
        --------------------------
        Args:
                * r: pd.DataFrame where index is geneSymbol and
                        rank_idx is the ranking index to use
                * gmt: path to gmt file (database to use for gsea)
                * rank_idx: column of dataframe to use for ranked enrichment
                * nperm: number of permutations for fgsea
                * nproc: how many processes to use for fgsea
                * min_size: minimum size of overlap for fgsea

        Returns:
                * fgsea_df: pd.DataFrame of results
        """
        fgsea = importr("fgsea")

        if rank_idx is None:
                rank_idx = r.columns[0]

        rrank = rpy2.robjects.vectors.FloatVector(r[rank_idx])
        rrank.names = list(r.index)

        fgsea_rresult = fgsea.fgsea(fgsea.gmtPathways(gmt), rrank, nperm=nperm, nproc=0, minSize=min_size)
        fgsea_df = pd.DataFrame.from_dict({key:np.asarray(fgsea_rresult.rx2(key)) for key in fgsea_rresult.names})

        try:
                fgsea_df['qval'],_ = qvalue(fgsea_df['pval'].values)
        except:
                try:
                        fgsea_df['qval'],_ = qvalue(fgsea_df['pval'].values, lambda_qvalue=0.5)
                except:
                        fgsea_df['qval'] = None

        return fgsea_df

# ----------------------------------
# Fisher Exact
# ----------------------------------
def build_counts_mat(df, cluster_id, cluster_label='consensus', description_label='grading'):
        """
        Build counts matrix.
        """
        df1 = pd.DataFrame(df[df[cluster_label] == cluster_id].groupby(description_label).size())
        df2 = pd.DataFrame(df[df[cluster_label] != cluster_id].groupby(description_label).size())
        full_df = pd.concat([df1,df2], 1)
        full_df.columns = ['in','out']

        full_df = full_df.fillna(0)
        full_df = full_df.astype(int)

        return full_df

def build_2x2(df):
        """Build 2x2 matrix."""
        df_out = df.sum(0) - df
        d = {}

        for i in df.index:
                cluster_i_df = pd.DataFrame(df.loc[i]).T
                cluster_o_df = pd.DataFrame(df_out.loc[i]).T

                cluster_i_df.index = ["i_cluster"]
                cluster_o_df.index = ["o_cluster"]

                d[i] = pd.concat((cluster_i_df, cluster_o_df))
        return d

def run_fisher_exacts(table_dict):
        """
        Run Fisher Exacts
        """
        from scipy.stats import fisher_exact

        indices = np.array(list(table_dict.keys()))
        odds_r = np.zeros(indices.shape[0])
        p_val = np.zeros(indices.shape[0])

        for i,idx in enumerate(indices):
                odds_r[i], p_val[i] = fisher_exact(table_dict[idx], alternative='greater')

        return pd.DataFrame(
                np.concatenate((odds_r[:,np.newaxis], p_val[:,np.newaxis]), axis=1),
                index=indices,
                columns=['odds_r','p_val']
        )

def compute_fisher_exact(
        labs,
        metadata_df,
        description_id="grading",
        label_id='consensus',
        fdr_alpha=0.05,
        fdr_method='fdr_bh'
        ):
        """
        Compute fisher exact.
        """
        from statsmodels.stats.multitest import multipletests

        fe_df = list()

        for lab in np.unique(labs[label_id]):
                lab_pval_df = run_fisher_exacts(
                        build_2x2(
                                build_counts_mat(metadata_df, lab, cluster_label=label_id, description_label=description_id)
                        )
                )
                lab_pval_df['id'] = lab
                fe_df.append(lab_pval_df)

        fe_df = pd.concat(fe_df).sort_values('p_val')
        _,fe_df['p_val_adj'],_,_ = multipletests(fe_df['p_val'], alpha=fdr_alpha, method=fdr_method)
        fe_df['qval'] = qvalue(fe_df['p_val'], lambda_qvalue=0.5)[0]

        return fe_df

def ptm_pval_fdr(res_df, protein_map, method='fdr_bh'):
    from statsmodels.stats.multitest import multipletests
    collapsed_res = pd.DataFrame(index=res_df['gene_name'].unique(), columns=res_df.columns)

    # For each protein, look at all PTM sites, select the one with the smallest p-value,
    # and multiply by the number of total PTM sites measured.
    for protein, protein_df in res_df.groupby('gene_name'):
        mostSig = protein_df['P.Value'].idxmin()
        collapsed_res.loc[protein] = res_df.loc[mostSig]
        collapsed_res.loc[protein,'P.Value'] *= protein_df.shape[0]
        if collapsed_res.loc[protein,'P.Value'] > 1:
            collapsed_res.loc[protein,'P.Value'] = 1
    # Run FDR on p-values
    collapsed_res['adj.P.Val'] = multipletests(collapsed_res['P.Value'], method=method, is_sorted=False,
                                              returnsorted=False)[1]

    # Storey q on p-values
    collapsed_res['qval'] = qvalue(collapsed_res['P.Value'])[0]
    # New GSEA rank on p-values
    collapsed_res['gsea_rank'] = collapsed_res.apply(lambda x: -np.log(x['adj.P.Val']) * x['logFC'], 1)
    return collapsed_res

# ----------------------------------
# Correlation
# ----------------------------------
def pull_geneset_from_gmt(gmt):
    """
    Get pathway genes from gene-set.
    """
    with open(gmt, 'r') as f:
        d = dict()
        for l in f.readlines():
            _l = l.strip().split("\t")
            d[_l[0]] = np.array(_l[2:])
    return d

from scipy.stats import spearmanr
def mcorr_global(x,y):
    """
    Subset for intersecting samples.
    """
    s = np.intersect1d(x.index, y.index)
    res = spearmanr(x[s], y[s])
    return pd.Series({'rho':res[0],'pval':res[1],'a':x.name,'b':y.name})

def m_corr_compute(a, b, genes, map_df, a_filt=None, b_filt=None, verbose=False):
    """
    Compute correlations w/ missimg values.
    -------------------------------
    Args:
        * a: assay_df (feats x samples)
        * b: assay_df (feats x samples)
        * genes: genes to include in analysis
        * map_df: dataframe mapping feats to geneSymbol
        * a_filt: min_frac of samples req. per feature
        * b_filt: min_frac of samples req. per feature
    """
    from tqdm import tqdm
    from multiprocessing import Pool, cpu_count

    def mcorr(x,y):
        """
        Subset for intersecting samples.
        """
        progBar.update()
        s = np.intersect1d(x.index, y.index)
        res = spearmanr(x[s], y[s])
        return pd.Series({'rho':res[0],'pval':res[1],'a':x.name,'b':y.name})

    class CorrTable:
        def __init__(self, a_df, b_df, a_prop, b_prop, progBar):
            self.a_df = a_df
            self.b_df = b_df
            self.a_prop = a_prop
            self.b_prop = b_prop
            self.res_idx = 0
            self.res_df = pd.DataFrame()
            self.pbar = progBar

        def postProcRes(self, map_df):
            self.res_df = self.res_df.set_index('a').join(map_df[['geneSymbol','feature']]).rename(
                columns={'geneSymbol':'a_gene','feature':'a_feature'}
                ).reset_index().rename(columns={'index':'a'})
            self.res_df = self.res_df.set_index('b').join(map_df[['geneSymbol','feature']]).rename(
                columns={'geneSymbol':'b_gene','feature':'b_feature'}
                ).reset_index().rename(columns={'index':'b'})
            self.res_df = self.res_df[['a','a_gene','a_feature','b','b_gene','b_feature','rho','pval','a_frac','b_frac']]
            return 0

        def runCorrMultiThread(self):
            def getCorr(self, a_feat):
                cpu_use = cpu_count()-1
                if cpu_use == 0:
                    cpu_use += 1
                print(f"Using {cpu_use} CPUs")
                p = Pool(cpu_use)
                r = p.starmap(mcorr_global, zip([self.a_df.loc[a_feat].dropna()]*self.b_df.shape[0],
                                                [self.b_df.loc[b_feat].dropna() for b_feat in self.b_df.index]))
                p.close()
                p.join()
                aFeat_df = pd.DataFrame(r)
                aFeat_df['a_frac'] = aFeat_df['a'].map(lambda x: self.a_prop.loc[x])
                aFeat_df['b_frac'] = aFeat_df['b'].map(lambda x: self.b_prop.loc[x])
                return aFeat_df
            for feature in self.a_df.index:
                self.res_df = pd.concat([self.res_df,getCorr(self, feature)])
                self.pbar.update()
            return 0

    feats = list(map_df[map_df['geneSymbol'].isin(genes)].index)
    feats += list(genes) # Includes Gene Centric Protein
    feats = np.unique(feats)

    a_idx = np.intersect1d(a.index, feats)
    b_idx = np.intersect1d(b.index, feats)

    if verbose: print(" * {} feats in A".format(a_idx.shape[0]))
    if verbose: print(" * {} feats in B".format(b_idx.shape[0]))

    _a = a.loc[a_idx].copy()
    _b = b.loc[b_idx].copy()
    _a_prop = _a.notna().sum(1)/_a.shape[1]
    _b_prop = _b.notna().sum(1)/_b.shape[1]

    if a_filt is not None:
        a_idx = _a_prop[_a_prop>a_filt].index
        if verbose: print("   * {} / {} kept (A)".format(a_idx.shape[0], _a.shape[0]))
        _a = _a.loc[a_idx]

    if b_filt is not None:
        b_idx = _b_prop[_b_prop>b_filt].index
        if verbose: print("   * {} / {} kept (B)".format(b_idx.shape[0], _b.shape[0]))
        _b = _b.loc[b_idx]

    progBar = tqdm(total=b_idx.shape[0], position=0, leave=True)
    corrTabObj = CorrTable(_a, _b, _a_prop, _b_prop, progBar)
    corrTabObj.runCorrMultiThread()
    corrTabObj.postProcRes(map_df)

    return corrTabObj.res_df

def corr_compute(df, type='spearman', drop_triangle=False):
    """
    Compute correlations with p-values for all columns in a matrix.
    Does not handle missing values.

        df: (samples x features)

    Adapted from: https://stackoverflow.com/questions/52741236/how-to-calculate-p-values-for-pairwise-correlation-of-columns-in-pandas
    """
    from scipy.stats import spearmanr, pearsonr

    def _collapse_tri(mat):
        """Collapse triangular matrix."""
        keep = np.triu(np.ones(mat.shape)).astype('bool').reshape(mat.size)
        if drop_triangle:
            return mat.stack()[keep].reset_index()
        else:
            return mat.stack().reset_index()

    mat = df.values.T
    K = len(df.columns)
    correl = np.empty((K,K), dtype=float)
    p_vals = np.empty((K,K), dtype=float)

    for i, ac in enumerate(mat):
        for j, bc in enumerate(mat):
            if i > j:
                continue
            else:
                if type=='spearman':
                    corr = spearmanr(ac,bc)
                elif type=='pearson':
                    corr = pearsonr(ac,bc)

            correl[i,j] = corr[0]
            correl[j,i] = corr[0]
            p_vals[i,j] = corr[1]
            p_vals[j,i] = corr[1]

    df_p = _collapse_tri(pd.DataFrame(p_vals, index=df.columns, columns=df.columns))
    df_corr = _collapse_tri(pd.DataFrame(correl, index=df.columns, columns=df.columns))

    df_p = df_p.rename(columns={0:'p-val'})
    df_corr = df_corr.rename(columns={0:'rho'})

    return pd.merge(df_p, df_corr).sort_values('rho')

# ----------------------------------
# Immune Deconvolution
# ----------------------------------
def loadImmuneDeconv(outdir: str):
    """
    Load all immune deconvolution results.
    ----------------------------
    Combines into dictionary.
    """
    # ESTIMATE
    estimate_df = pd.read_csv(os.path.join(outdir, "estimate", "cptac_tpm_score.gct"), skiprows=2, sep="\t", index_col=0).iloc[:,1:].T

    # xCell Transcriptome
    xcell_rna_df = pd.read_csv(os.path.join(outdir, "xCell", "rna_xcell.tsv"), sep='\t', index_col=0).T
    xcell_rna_scores_df = xcell_rna_df[['ImmuneScore','StromaScore', 'MicroenvironmentScore']]
    xcell_rna_df = xcell_rna_df.drop(columns=['ImmuneScore','StromaScore', 'MicroenvironmentScore'])

    # xCell Proteome
    xcell_prot_df = pd.read_csv(os.path.join(outdir, "xCell", "prot_xcell.tsv"), sep='\t', index_col=0).T
    xcell_prot_scores_df = xcell_prot_df[['ImmuneScore','StromaScore', 'MicroenvironmentScore']]
    xcell_prot_df = xcell_prot_df.drop(columns=['ImmuneScore','StromaScore', 'MicroenvironmentScore'])

    # CIBERSORT
    ciber_df = pd.read_csv(os.path.join(outdir,"cibersort", "cibersort.tsv"), sep='\t').iloc[:,:-3]

    # Immune Subtype Classifier
    isc_names = {
        "1": "Wound Healing",
        "2": "IFN-Gamma",
        "3": "Inflammatory",
        "4": "Lymphocyte Depleted",
        "5": "Immunologically Quiet",
        "6": "TGF-Beta"
    }

    immune_sc_df = pd.read_csv(os.path.join(
        outdir,
        "immune_subtype_classifier",
        "isc_results.tsv"
    ), sep='\t', index_col=0).set_index("SampleIDs")

    immune_sc_df = immune_sc_df.rename(columns=isc_names)
    immune_sc_df['BestCall'] = immune_sc_df['BestCall'].apply(lambda x: isc_names[str(x)])

    # Combine
    s_o = xcell_rna_df.index

    resd = dict()
    resd['xcell_rna'] = xcell_rna_df
    resd['xcell_rna_scores'] = xcell_rna_scores_df.loc[s_o,:]
    resd['xcell_prot'] = xcell_prot_df.loc[s_o,:]
    resd['xcell_prot_scores'] = xcell_prot_scores_df.loc[s_o,:]
    resd['estimate'] = estimate_df.loc[s_o,:]
    resd['ciber'] = ciber_df.loc[s_o,:]
    resd['isc'] = immune_sc_df.loc[s_o,:]

    return resd

def load_ssgsea_results(rna_file, prot_file):
    """
    Load ssGSEA results.
    """
    ssgsea_rna = gen_ptmgsea_df(GCT(rna_file, parse_fn=False))
    ssgsea_prot = gen_ptmgsea_df(GCT(prot_file, parse_fn=False))

    ssgsea_rna['feature'] = 'transcriptome'
    ssgsea_prot['feature'] = 'proteome'

    # Full file
    ssgsea_df = pd.concat([ssgsea_rna,ssgsea_prot])
    ssgsea_df = ssgsea_df.set_index("id")
    ssgsea_df["NES"] = ssgsea_df["NES"].astype(float)
    ssgsea_df["fdr.pvalue"] = ssgsea_df["fdr.pvalue"].astype(float)

    # Matrices
    ssgsea_rna_mat = ssgsea_df[ssgsea_df['feature']=='transcriptome'].reset_index().pivot(index=['id'], columns=['pathway'], values=['NES'])
    ssgsea_prot_mat = ssgsea_df[ssgsea_df['feature']=='proteome'].reset_index().pivot(index=['id'], columns=['pathway'], values=['NES'])
    ssgsea_rna_mat.columns = ssgsea_rna_mat.columns.droplevel()
    ssgsea_prot_mat.columns = ssgsea_prot_mat.columns.droplevel()

    return ssgsea_df, ssgsea_rna_mat, ssgsea_prot_mat
