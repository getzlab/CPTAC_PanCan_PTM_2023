import pandas as pd
import sys
import argparse
sys.path.append("../funcs/")
from proteomics import GCT

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

def main():
    parser = argparse.ArgumentParser(description='Save & write out ptm-gsea ssGSEA summary.')
    parser.add_argument('-i', '--input', help='<Required> Input path to PTM-GSEA combined output.')
    parser.add_argument('-o', '--output', help='<Required> Output file.')
    parser.add_argument('-f', '--formatstr', default='', help='Format naming str.')

    args = parser.parse_args()

    res_df = gen_ptmgsea_df(GCT(args.input, parse_fn=False), fmt=args.formatstr)
    res_df.to_csv(args.output, sep='\t')

if __name__ == "__main__":
	main()
