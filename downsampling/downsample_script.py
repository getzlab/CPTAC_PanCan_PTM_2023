import pandas as pd
import signatureanalyzer as sa
from typing import Union
import os

def run_downsampling_ardnmf(
    X: pd.DataFrame,
    n_samples: int,
    n_iter: int = 1,
    seed: Union[None, int] = 42,
    out_dir: str = "."
    ):
    """
    Run downsampling analysis.

    Parameters:
        X: input for ARD-NMF.
        n_samples: number to sample without replacement
        seed: random seed
    """
    import numpy as np

    if seed is not None:
        np.random.seed(seed)

    samples = np.array(X.columns)

    # Output directory
    for n in range(n_iter):
        idx = np.random.choice(samples, n_samples, replace=False)

        sa.run_matrix(
            X.loc[:,idx],
            outdir=os.path.join(out_dir, "sample_{}".format(n)),
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
output_dir = "downsamples"
n_samples = [100,200,300,400,500,600,700,800,900,1000]
n_iter = 30

# Run for downsamples
for n in n_samples:
    run_downsampling_ardnmf(X, n, n_iter=n_iter, out_dir = os.path.join(output_dir, "n_{}".format(n)))
