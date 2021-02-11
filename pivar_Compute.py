from __future__ import division
import string,re,sys,copy, math,collections, itertools


def ProcessCLI(args):
    
    contiglen=readidx(args[2])
    getPiVar(args[1],contiglen,args[3])

def getN2(n):
    return n*(n-1)/2

def getN3(n):
    return n*(n-1)*(n-2)/2

def CalcSitehsVar(snpinfo):
    
    pos=snpinfo[1]
    n=float(snpinfo[4])
    maf=float(snpinfo[-1])
    n2=getN2(n);n3=getN3(n)
    hs=2*maf*(1-maf)
    Vars= hs*(1-hs +(1-2*hs)*n3/n2)/n2
    
    return pos,hs,Vars

def getPiVar(vcflike,contiglen,out):

    sitehet=open(out+".hs.out","w")
    ContigHet= open(out+".pivar.out","w")

    contigfreq={}
    
    for line in open(vcflike):
        
        words=line.strip().split(None)
        pos,hs,var=CalcSitehsVar(words)
        
        if words[0] not in contigfreq.keys():
            contigfreq[words[0]]=[]
            contigfreq[words[0]]=[(hs,var)]
        else:
            contigfreq[words[0]]+= [(hs,var)]
            
        sitehet.write(words[0] + "\t" + str(pos) + "\t" + "{0:.4f}".format(hs) + "\t" + "{0:.4f}".format(var) + "\n")

    allcontigs=[]
    
    for k in contigfreq.keys():
        
        allloci=[]
        for pos in contigfreq[k]:
            allloci+=[pos] 
            
        if len(allloci)>0:
            allcontigs+=[(k,[sum(x) for x in zip(*allloci)])]

    for c in allcontigs:
        
        ContigHet.write(c[0] + "\t"+ str(contiglen[c[0]])+"\t" + "{0:.4f}".format(c[1][0]) + "\t" + "{0:.4f}".format(c[1][1]) + "\t" +"{0:.4f}".format(c[1][0]/contiglen[c[0]])+ "\n")

def readidx(idx):

    contiglen={}
    
    for line in open(idx,"r"):
        
        words=line.strip().split(None)
        try:
            contiglen[words[0]]=int(words[1])
        except KeyError :
            pass
        
    return contiglen 

if __name__ == "__main__":
    ProcessCLI(sys.argv)

