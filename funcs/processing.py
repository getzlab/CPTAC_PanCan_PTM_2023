import pandas as pd
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import glob
import seaborn as sns
from scipy.stats import norm

import warnings
warnings.filterwarnings('ignore')

from plotting import CPTAC_CMAP

def print_d(x):
    """Print func"""
    for k in x.keys():
        print("{}: {} samples x {} features.".format(k, x[k].shape[0], x[k].shape[1]))

def cohort_aggr(file_paths: list, sites=None):
    """
    Aggregate load proteomics files into single file.
    """
    X_list = list()
    meta_list = list()

    for file in file_paths:
        X, meta = pload(file)
        meta.loc['cohort'] = file.split("/")[-1].split("_")[0].upper()
        meta_list.append(meta)

        if sites is not None:
            sites_intersect = set.intersection(set(X.index),sites)
            missing_sites = set(sites)-set(sites_intersect)
            X = pd.concat([X,pd.DataFrame(None,columns=X.columns, index=missing_sites)])
            X_list.append(X.loc[sites])
        else:
            X_list.append(X)

    X = pd.concat(X_list,1).astype(float).T
    meta_df = pd.concat(meta_list,1)

    return X, meta_df.T

def plot_scatter(x, ax, axis=0, c=None, cmap=None):
    """
    Plot scatterplot.
    ----------------------
    Args:
        * x: input dataframe
        * ax: axis to plot on
        * axis: axis to compute std. & mean over
        * c: color (pd.Series)
    """
    if c is not None:
        cohorts = list(set(c))
        cohorts = np.sort(cohorts)

        if cmap is not None:
            color_dict = {cohorts[i]:x for i,x, in enumerate(sns.color_palette(cmap, len(set(c))))}
        else:
            color_dict = CPTAC_CMAP

        ax.scatter(x.mean(axis), x.std(axis), alpha=0.4, c=c.apply(lambda x: color_dict[x]))

        xmin, xmax = ax.get_xlim()
        ymin, ymax = ax.get_ylim()

        for k,v in CPTAC_CMAP.items():
            if k in cohorts:
                m = ax.scatter(xmin-1, ymin-1, alpha=.8, c=np.array(v)[np.newaxis,:], label=k)

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.legend()

    ax.scatter(x.mean(axis), x.std(axis), alpha=0.4)
    ax.set_xlabel("Mean", fontsize=14)
    ax.set_ylabel("Std. Dev.", fontsize=14)

def plot_mean_std_grid(X, meta, figsize=(12,8), y1_lim=(0,12), y2_lim=(0,5.5)):
    """
    Plots grid of mean & std dev plots for each cohort.
    -----------------------
    Args:
        * X: dict of matrices
        * meta: dict of metadata
        * figsize: size of figure
        * y1_lim: y_lim for site std. dev.
        * y2_lim: y_lim for patient std. dev.

    Returns:
        None
    """
    fig,axes = plt.subplots(2,len(X.keys()), figsize=figsize, sharey=False)

    for idx,feat in enumerate(X.keys()):
        plot_scatter(X[feat], axes[0,idx])
        plot_scatter(X[feat], axes[1,idx], axis=1, c=meta[feat]['cohort'])
        axes[0,idx].set_title(feat.capitalize(), fontsize=16)
        axes[0,idx].set_xlabel(None)

        if idx > 0:
            axes[0,idx].set_ylabel(None)
            axes[1,idx].set_ylabel(None)

        axes[1,idx].legend(frameon=False)

        axes[0,idx].set_ylim(y1_lim)
        axes[1,idx].set_ylim(y2_lim)

    axes[0,0].set_ylabel("Std. Dev. by Site", fontsize=14)
    axes[1,0].set_ylabel("Std. Dev. by Patient", fontsize=14)

def filter_NA(X, thresh=0.3, axis=0, ax=False, title="", meta=None, cmap=None):
    """
    Filter out NA values below certain threshold.
    ------------------------------
    Args:
        * X: pd. DataFrame input
        * thresh: percent threshold across that axis
        * axis: axis to use for threshold
        * ax: axis to plot; if not provided, plot will not be created
        * meta: pd.Series of cohort information to split by
            * NOTE: to select what sites to keep, we will take the union
                    across all sites selected for each cohort.

    Returns:
        * np.ndarray sites to keep
    """
    if meta is not None:
        sites_keep_by_cohort = {}
        na_proportions_by_cohort = {}

        for m in np.unique(meta):
            na_proportions = X.loc[meta[meta==m].index].isnull().sum(axis=axis) / X.loc[meta[meta==m].index].shape[axis]
            sites_keep = na_proportions[na_proportions < thresh]
            sites_keep = sites_keep.index.values
            np.sort(sites_keep)
            sites_keep_by_cohort[m] = sites_keep
            na_proportions_by_cohort[m] = na_proportions

        # Return union of all kept sites across cohorts
        # sites_keep = np.array(list(set.intersection(*[set(x) for x in sites_keep_by_cohort.values()])))
        sites_keep = np.array(list(set.union(*[set(x) for x in sites_keep_by_cohort.values()])))
        np.sort(sites_keep)

        # -----------------
        # Plot
        # -----------------
        if ax:
            cohorts = np.unique(meta)
            if cmap is None:
                color_dict = CPTAC_CMAP
            else:
                color_dict = {cohorts[i]:x for i,x, in enumerate(sns.color_palette(cmap, cohorts.shape[0]))}

            for m in np.unique(meta):
                ax.hist(na_proportions_by_cohort[m], edgecolor='black', bins=25, color=color_dict[m], alpha=0.2)
                #sns.distplot(na_proportions_by_cohort[m], ax=ax, kde_kws={'bw':0.05})

            ax.set_xlim(0,1)
            ax.set_ylabel("Freq", fontsize=(14))
            ax.set_title("{}{} / {} Sites".format(title, sites_keep.shape[0], na_proportions.shape[0]), fontsize=16)

            ax.legend(np.unique(meta), frameon=False)
            ax.axvline(thresh, color='black', linestyle='-.')

    else:
        na_proportions = X.isnull().sum(axis=axis) / X.shape[axis]
        sites_keep = na_proportions[na_proportions < thresh]
        sites_keep = sites_keep.index.values
        np.sort(sites_keep)

        # -----------------
        # Plot
        # -----------------
        if ax and meta is None:
            ax.hist(na_proportions, edgecolor='black', color='lightblue', bins=int(na_proportions.shape[0]/200))

            ax.set_xlim(0,1)
            ax.set_ylabel("Freq", fontsize=(14))
            ax.set_title("{}{} / {} Sites".format(title, sites_keep.shape[0], na_proportions.shape[0]), fontsize=16)

            ax.axvline(thresh, color='black', linestyle='-.')
            ax.legend(['Filter Threshold'])

    return sites_keep

def clean_tumor_normal_annot(row):
    """
    This function annotates all samples as either tumor samples
    or normal samples from the varied proteomis metadata.
    """
    if row["Type"] in ("Normal", "NAT", "Adjacent_Normal"):
        return "Normal"
    elif "IR" in row.name:
        return "IR"
    elif row["Type"] == "Tumor":
        return "Tumor"
    else:
        if "Normal" in row.name:
            return "Normal"
        else:
            return "Tumor"

    return row["Type"]

def aggr_feature_cohorts(d: dict, intersect_patients: bool = True):
    """
    Takes a dictionary of matrices & intersects them to great
    one aggregated, cross-feature cohort.
    -------------
    Args:
        * d: dictionary mapping feature type --> pd.DataFrame (samples x features)
        * intersect_patients: bool --> use only intersecting patients

    Returns:
        * pd.DataFrame: (samples x features) aggregated
    """
    feature_types = np.sort(list(d.keys()))
    df_list = [d[ft].T for ft in np.sort(list(d.keys()))]

    if intersect_patients:
        return pd.concat(df_list,0).loc[:,set.intersection(*map(set, df_list))].drop_duplicates().T
    else:
        return pd.concat(df_list,0).T
