import pandas as pd
import signatureanalyzer as sa
from typing import Union
import os
import glob
import numpy as np
from tqdm import tqdm

def run_cohorts_ardnmf(
    X: pd.DataFrame,
    cohorts: list,
    clust_df: pd.DataFrame,
    seed: Union[None, int] = 42,
    out_dir: str = "."
    ):
    """
    Run downsampling analysis.

    Parameters:
        X: input for ARD-NMF.
        cohorts: list of groupings
        seed: random seed
    """
    import numpy as np
    import itertools

    if seed is not None:
        np.random.seed(seed)

    combs = list()

    for x in range(2,len(cohorts)):
        combs += list(itertools.combinations(cohorts,x))

    ######### Completed Runs #########
    _combs = list()

    for comb in combs:
        if not os.path.isdir(os.path.join(out_dir, '_'.join(np.sort(comb)))):
            _combs.append(comb)

    if len(_combs) > 0:
        print("{} / {} runs found - resuming remainder {}.".format(
            len(combs)-len(_combs), len(combs), len(_combs)))
        combs = _combs
    ##################################

    # Output directory
    for comb in tqdm(combs):
        idx = np.array(clust_df[clust_df['cohort'].isin(comb)].index)
        sa.run_matrix(
            X.loc[:,idx],
            outdir=os.path.join(out_dir, '_'.join(np.sort(comb))),
            nruns = 20,
            verbose=False,
            plot_results=False,
            K0=50,
            max_iter=10000,
            prior_on_W='L2',
            prior_on_H='L2',
            objective='gaussian'
        )

# Load Data
X = pd.read_parquet("../../data/processed/061721/imputed_res/pan_reg_X.parquet")
output_dir = "downsamples_by_cohort"

os.makedirs(output_dir, exist_ok=True)

clust_df = pd.read_csv('../signatures/061721_imputed_res_reg/mappings_pancan_102322.tsv', sep='\t',index_col=0)
cohorts = list(set(clust_df['cohort']))

# Run for downsamples
run_cohorts_ardnmf(X, cohorts, clust_df, out_dir = output_dir)
