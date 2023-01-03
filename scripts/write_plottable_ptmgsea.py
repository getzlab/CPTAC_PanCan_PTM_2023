import pandas as pd
import sys
import argparse

def ptmgsea_plottable(setFeature, df, setIDTrue=None, setIDFalse=None):
    df['feature'] = setFeature
    if setIDTrue and setIDFalse:
        df['id'] = df['id'].map(lambda x: setIDTrue if x==True else setIDFalse)
    df = df.rename(columns={'fdr.pvalue':'padj'})
    return df

def main():
    parser = argparse.ArgumentParser(description="Make PTM-GSEA output plottable with existing code.")
    parser.add_argument('--raw', required=True, help='Input path to PTM-GSEA raw phosphoproteome combined output.')
    parser.add_argument('--res', required=True, help='Input path to PTM-GSEA residualized phosphoproteome combined output.')
    parser.add_argument('--trueID', required=False, default=None, help='ID to set for True')
    parser.add_argument('--falseID', required=False, default=None, help='ID to set for False')
    parser.add_argument('-o', '--outpath', required=True, help='Output file')
    
    args = parser.parse_args()

    raw_df = pd.read_csv(args.raw, sep='\t', index_col=0)
    res_df = pd.read_csv(args.res, sep='\t', index_col=0)
    
    pd.concat([ptmgsea_plottable('Phosphoproteome Raw', raw_df, args.trueID, args.falseID),
               ptmgsea_plottable('Phosphoproteome Res.', res_df, args.trueID, args.falseID)]).to_csv(args.outpath, sep='\t')

if __name__ == '__main__':
    main()
    
