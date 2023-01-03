import pandas as pd
import subprocess
from subprocess import Popen, PIPE
import os
import numpy as np
import re

def _run(x):
    "Generic run."
    res = subprocess.run(x.split(' '), stdout=PIPE, stderr=PIPE, universal_newlines=True)
    res.stdout = res.stdout.strip().split('\n')
    res.stderr = res.stderr.strip().split('\n')
    return res

def load_de(x):
    """Load de result."""
    df = pd.read_csv(x, sep='\t', index_col=0).reset_index()
    df = df[df['id']]
    df.loc[:,'id'] = x.split("/")[-2].split("_")[1]
    return df

def load_immune_de(x):
    """Load de result."""
    df = pd.read_csv(x, sep='\t', index_col=0).reset_index()
    if 'Score' in x:
        df = df[df['id']]
        df.loc[:,'id'] = x.split("/")[-2].split("_")[1]
    return df

def load_mut_de(x):
    """Load de result."""
    df = pd.read_csv(x, sep='\t', index_col=0).reset_index()
    df = df[df['id']]
    df.loc[:,'id'] = x.split("/")[-2]
    return df

def write_log(fout,res):
    """Write output logs."""
    with open(fout+'.out', "w") as f:
        for l in res.stdout:
            f.write(l+'\n')

    with open(fout+'.err', "w") as f:
        for l in res.stderr:
            f.write(l+'\n')

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

def gen_phospho_ptmgsea(diffexp_file, map_df, outfile, weights='gsea_rank', feature='phosphoproteome'):
    """
    Generate input file for run name.
    """
    from ast import literal_eval
    
    if isinstance(diffexp_file, str):
        de_df = pd.read_csv(diffexp_file, sep='\t', index_col=0)
    else:
        de_df = diffexp_file

    de_df = de_df[de_df['feature']==feature]

    de_df = de_df.loc[:,['gene_name', 'logFC', 'AveExpr', 't', 'P.Value', 'adj.P.Val', 'B','qval', 'id', weights]]
    de_df = de_df.join(map_df)
    de_df = de_df[~de_df['VMsiteFlanks'].isna()]

    def selectRank(df, weights):
        # input dataframe: stats for a given gene_site_flank
        # Return: pd.Series for the peptide to use
        if df['single_ptm'].any():
            max_idx = df[df['single_ptm']][weights].abs().idxmax()
            return df.loc[max_idx]
        else:
            max_idx = df[weights].abs().idxmax()
            return df.loc[max_idx]
    
    phosph_var_agg = pd.DataFrame()
    phosph_X = pd.DataFrame()
    first = True
    for condition, cond_df in de_df.groupby('id'):
        phosph_var = cond_df[['geneSymbol','id.description','accession_number','protein_mw','variableSites','sequence',
                              'sequenceVML','VMsiteFlanks','causalpath_adjusted_id', weights]]
        phosph_var['VMsiteFlanks_l'] = phosph_var['VMsiteFlanks'].apply(lambda x: x.replace('|', "', '")).apply(literal_eval)
        phosph_var['varSites_l'] = phosph_var['causalpath_adjusted_id'].apply(lambda x: re.findall(r'[STY][0-9]*[sty]', x))
        phosph_var = phosph_var[phosph_var['VMsiteFlanks_l'].map(len) == phosph_var['varSites_l'].map(len)]
        phosph_var['gene_uniprotPos_flank_l'] = phosph_var.apply(lambda x: [x['geneSymbol'] + '_' + x['varSites_l'][i] + 
                                                                            '_' + x['VMsiteFlanks_l'][i] 
                                                                            for i in range(len(x['varSites_l']))], 1)
        phosph_var = phosph_var.explode('gene_uniprotPos_flank_l')
        phosph_var['single_ptm'] = phosph_var['varSites_l'].map(len) == 1

        phosph_var = phosph_var.reset_index().groupby('gene_uniprotPos_flank_l').apply(lambda x: selectRank(x, weights))
        phosph_X_temp = pd.DataFrame(phosph_var[weights])
        phosph_X_temp.columns = [condition]
        phosph_X = pd.concat([phosph_X, phosph_X_temp], 1)
        
        var_mergeOn = ['gene_uniprotPos_flank_l','geneSymbol']
        var_metaKeep = ['gene_uniprotPos_flank_l','geneSymbol','index']
        phosph_var = phosph_var.drop(columns=['gene_uniprotPos_flank_l']).reset_index()
        phosph_var = phosph_var[var_metaKeep]
        var_suffix_cols = [f"{x}_{condition}" for x in phosph_var.columns if x not in var_mergeOn]
        phosph_var.columns = var_mergeOn + var_suffix_cols

        if first:
            phosph_var_agg = phosph_var
            first = False
        else:
            phosph_var_agg = pd.merge(phosph_var_agg,
                                      phosph_var,
                                      on=var_mergeOn,
                                      how='outer')

    phosph_var_agg['ptmGSEA'] = phosph_var_agg['gene_uniprotPos_flank_l'].map(lambda x: x.split('_')[-1].replace('-','_').upper() + '-p')
    phosph_var_agg = phosph_var_agg.set_index("ptmGSEA")
    phosph_X.index = phosph_X.index.map(lambda x: x.split('_')[-1].replace('-', '_').upper() + '-p')
    phosph_X.index.name = 'ptmGSEA'
    
    phosph_obs = pd.DataFrame(phosph_X.columns, columns=['S']).set_index("S")
    
    if isinstance(diffexp_file, str):
        phosph_obs['run_name'] = diffexp_file

    # Write GCT
    write_gct(phosph_X, phosph_obs, phosph_var_agg, outfile)

def run_differential_expression(formats, scripts_dir, verbose=False):
    """
    Run differential expression.
    """
    # ----------- Run Diffexp
    cmd = "Rscript {}/run_limma.R".format(scripts_dir)

    if 'covar' in formats:
        cmd += " --covar {}".format(formats['covar'])
    if 'proteome' in formats:
        cmd += " --proteome {}".format(formats['proteome'])
    if 'phosphoproteome' in formats:
        cmd += " --phosphoproteome {}".format(formats['phosphoproteome'])
    if 'phosphoproteome_res' in formats:
        cmd += " --phosphoproteome_res {}".format(formats['phosphoproteome_res'])
    if 'transcriptome' in formats:
        cmd += " --transcriptome {}".format(formats['transcriptome'])
    if 'acetylome' in formats:
        cmd += " --acetylome {}".format(formats['acetylome'])
    if 'acetylome_res' in formats:
        cmd += " --acetylome_res {}".format(formats['acetylome_res'])
    if 'subset' in formats:
        cmd += ' --subset {}'.format(formats['subset'])
    if 'covar_to_use' in formats:
        cmd += ' --covar_to_use {}'.format(formats['covar_to_use'])
    if 'minObs' in formats:
        cmd += ' --minObs {}'.format(formats['minObs'])
    cmd += ' --labels {}'.format(formats['input']) \
        +' --label_id {}'.format(formats['clust']) \
        +' --feature_maps {}'.format(formats['prot_maps']) \
        +' --output {}'.format(formats['out_limma'])

    os.makedirs(formats['out_limma'], exist_ok=True)

    if verbose: print(cmd)
    res = _run(cmd)
    write_log(os.path.join(formats['out_limma'], "de.log"), res)

    # ----------- Post Process Diffexp
    cmd = "python3 {}/postprocess_limma_de.py".format(scripts_dir)
    cmd += ' -i {}'.format(formats['out_limma'])\
        +' -f {}'.format(formats['prot_maps'])\
        +' -o {}'.format(formats['out_limma'])

    if verbose: print(cmd)
    res = _run(cmd)
    write_log(os.path.join(formats['out_limma'], "de.postprocess.log"), res)

    # ----------- Prepare PTM GSEA
    map_df = pd.read_csv(formats['prot_maps'], sep='\t', index_col=0)
    gen_phospho_ptmgsea(
        os.path.join(formats['out_limma'], "full_diffexp_results.tsv"),
        map_df,
        os.path.join(formats['out_limma'], "de_phosph_raw_ptmgsea_p.gct"),
        weights='gsea_rank_p',
        feature='phosphoproteome'
    )
    gen_phospho_ptmgsea(
        os.path.join(formats['out_limma'], "full_diffexp_results.tsv"),
        map_df,
        os.path.join(formats['out_limma'], "de_phosph_res_ptmgsea_p.gct"),
        weights='gsea_rank_p',
        feature='phosphoproteome_res'
    )
