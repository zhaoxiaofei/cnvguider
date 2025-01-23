#!/usr/bin/env python

import sys
import pandas as pd

def weighted_median(df, val, weight):
    df_sorted = df.sort_values(val)
    cumsum = df_sorted[weight].cumsum()
    cutoff = df_sorted[weight].sum() / 2.
    return df_sorted[cumsum >= cutoff][val].iloc[0]

def main():
    bedtable1 = pd.read_csv(sys.stdin, sep='\t', header=0)
    bedtable1['region_size'] = bedtable1['end_37'] - bedtable1['start_37']
    medianDP = weighted_median(bedtable1, 'obsDP', 'region_size')
    #print(medianDP)
    bedtable2 = bedtable1.copy()
    bedtable2['obsCN'] = (bedtable1['obsDP'] / medianDP).round(0).astype(int)
    bedtable2[['#chr_37', 'start_37', 'end_37', 'obsCN']].to_csv(sys.stdout, sep='\t', header=True, index=None)

if __name__ == '__main__': main()

