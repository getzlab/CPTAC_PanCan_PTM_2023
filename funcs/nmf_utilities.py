import sys
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
from typing import Union
from sys import stdout
import signatureanalyzer as sa

def nmf_loader(filepath: str):
    """
    NMF loading function.
    ----------------------------
    Args:
        * path to dir nmf_output.h5 file

    Returns:
        * tuple (pd. DataFrame)
            X, H, W, signatures, markers
    """
    X = pd.read_hdf(filepath,'X')
    X.loc[[x for x in X.index if x.endswith('_n')]] = X.loc[[x for x in X.index if x.endswith('_n')]]*-1

    H = pd.read_hdf(filepath,'H')
    W = pd.read_hdf(filepath,'W')
    signatures = pd.read_hdf(filepath,'signatures')
    markers = pd.read_hdf(filepath, 'markers')

    return X, H, W, signatures, markers

def pivot_props(df: pd.DataFrame, x='cohort', y='clusters', norm=True):
    """
    Pivot Proportions.
    ------------------
    Args:
        * df: pd.DataFrame input
        * x: x-value of combined dataframe
        * y: y-value of combined dataframe
        * norm: whether or not to convert to proportion

    Return:
        * combined and pivoted dataframe
            ex. generates proportions of cohorts over each cluster id
    """
    assert x in df and y in df, "Ensure x & y are columns in input dataframe."

    pivot_df = pd.DataFrame(df[[x,y]].groupby([x,y]).size()).reset_index().rename(columns={0:'count'}).pivot(
        index=y, columns=x, values='count'
    ).fillna(0).astype(int)

    if norm:
        return pivot_df.T.div(pivot_df.sum(1)).T
    else:
        return pivot_df
