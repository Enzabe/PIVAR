
import pandas as pd
import numpy as np 
import sys

def ProcessCLI(args):
     
    df=getdf(args[1])
    df.columns=["contig","Len","MapReads","UMapReads"]
    df=df[df.contig !='*']
    df=df[df.MapReads >1000]    
    #df=df.set_index(['contig','Len'])
    #df2=df.copy()
    #df=df.T
    #df=df.div(df.sum(axis=1), axis=0)
    #df=df.T
    #df=df.reset_index()
    df['RPKM']=df['MapReads']/df['MapReads'].sum()
    df['RPKM'] = df.apply(lambda row: row.RPKM*1000000/row.Len, axis=1)


    df=df.drop(['UMapReads'],axis=1)
    df=df.set_index(['contig','Len'])
    df=df.round(4)
    #df2['RPKM']=df['RPKM']
    #df2=df2.drop(['UMapReads'],axis=1)
    #df=df2
    #print(df)
    #return
    df.to_csv(args[2],sep=' ')
    
def getdf(indf):
    df=pd.read_csv(indf,sep="\t", header=None)
    return df 


if __name__ == "__main__":
    ProcessCLI(sys.argv)
