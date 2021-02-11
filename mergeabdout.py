import pandas as pd
import sys
import numpy as np
from collections import defaultdict
from collections import Counter


def ProcessCLI(args):
    dfs=getfile(args[1])
    mergedf=mergedfs(dfs)
    mergedf=mergedf.round(4)
    mergedf.to_csv(args[2],sep=' ',na_rep='NA')
    
def getfile(infile):
    afiles=[ line.strip()+'.abd' for line in open(infile,'r').readlines()]
    return afiles 
       
def mergedfs(abdfiles):
    
    df=getdf(abdfiles[0])
    if len(abdfiles)>1:
        for i,k in enumerate(abdfiles[1:]):
            df=df.merge(getdf(k),left_index=True, right_index=True,how='outer',suffixes=('', '_'+str(i+1)))
    return df 


def getdf(indf):
    df=pd.read_csv(indf,sep=" ", header=0)
    df=df.set_index(['contig','Len'])
    df=df.drop(columns=['MapReads'])
    return df 

if __name__ == "__main__":
    ProcessCLI(sys.argv)
