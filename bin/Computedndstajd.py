from __future__ import division
import string,re,sys,copy,collections
import ctypes 
from math import * 
import numpy as np  

def ProcessCLI(args):

    seqs=readseq(args[1])
    gids2=getgids(args[2])
    out=open(args[3],"w")

    codonns={}
    codontable={
        'ATT':'I', 'ATC':'I', 'ATA':'I',
        'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'TTA':'L', 'TTG':'L',
        'GTT':'V', 'GTC':'V', 'GTA':'V','GTG':'V',
        'TTT':'F', 'TTC':'F',
        'ATG':'M',
        'TGT':'C', 'TGC':'C',
        'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
        'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
        'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
        'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
        'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
        'TAT':'Y', 'TAC':'Y',
        'TGG':'W',
        'CAA':'Q', 'CAG':'Q',
        'AAT':'N', 'AAC':'N',
        'CAT':'H', 'CAC':'H',
        'GAA':'E', 'GAG':'E',
        'GAT':'D', 'GAC':'D',
        'AAA':'K', 'AAG':'K',
        'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
        'TAA':'STOP', 'TAG':'STOP', 'TGA':'STOP'
        }
    for codon in codontable.keys():
        n,s=getndns(codon,codontable)
        codonns[codon]=[n,s]
    
    gids={}
    
    for k in gids2.keys():
        if len(gids2[k].keys())>1:
            gids[k]=gids2[k]

    for gene in gids.keys():
        
        seq=seqs[">"+gene]
        dndstjd=calcdNdStajd(seq,gids,codontable,gene)
        gdt='\t'.join([gene]+[str(i) for i in dndstjd])
        out.write(gdt+"\n")

def calcdNdStajd(seq,gids,codontable,gene):
    
    N,S=calulateNandS(seq,codontable)
    pn=None 
    ps=None
    tajd=None
 
    if gene[-1]=="+":
        pn,ps,tajd=getpnpsfwd(seq,gids,codontable,gene)
        
    if gene[-1]=="-":
        pn,ps,tajd=getpnpsrevstr(seq,gids,codontable,gene)
    
    dN=(-3.0/4)*np.log(1-4*pn/3)
    dS=(-3.0/4)*np.log(1-4*ps/3)
    if dS>0:
        return [ round(N,4)+0.0,round(S,4)+0.0,round(pn,4)+0.0,round(ps,4)+0.0,round(dN,4)+0.0,round(dS,4)+0.0,round(dN/dS,4)+0.0,round(tajd,4)+0.0]
    else:
        return [ round(N,4)+0.0,round(S,4)+0.0,round(pn,4)+0.0,round(ps,4)+0.0,round(dN,4)+0.0,round(dS,4)+0.0,'na',round(tajd,4)+0.0]


def getpnpsrevstr(seq,gids,codontable,gene):
    ns=0
    nd=0
    r= len(translateDNA(seq,codontable))
    freqs=[]
    n=[]
    sites=len(gids[gene].keys())
    
    if gene[-1]=="-":
        seq=getReverseStrand(seq)

    if len(gids[gene].keys())>1: 
        
        for pos in gids[gene].keys():
            freqs.append(float(gids[gene][pos][5]))
            n.append(float(gids[gene][pos][2]))

            p=len(seq)-(int(pos)+1)
   
            if (p+1)%3==0:
                 
                try:
                    codon=seq[p-2:p+1]
                    ni,si=getndns(codon,codontable)
                    mutcodon= substbase(codon,2,getReverseStrand(gids[gene][pos][1]))
                    if codontable[codon]==codontable[mutcodon]:
                        ns+=(float(gids[gene][pos][4])/si)/float(gids[gene][pos][2])
                    if codontable[codon]!=codontable[mutcodon]:
                        nd+=(float(gids[gene][pos][4])/ni)/float(gids[gene][pos][2])
                except ZeroDivisionError as detail:
                    pass
                    #print('Error: on this ',gene,pos,detail)
                    
            if (p+1)%3==1:

                try:
                    
                    codon=seq[p:p+3]
                    ni,si=getndns(codon,codontable)
                    mutcodon= substbase(codon,0,getReverseStrand(gids[gene][pos][1]))

                    if codontable[codon]==codontable[mutcodon]:

                        ns+=(float(gids[gene][pos][4])/si)/float(gids[gene][pos][2])
                    elif codontable[codon]!=codontable[mutcodon]:
                        nd+=(float(gids[gene][pos][4])/ni)/float(gids[gene][pos][2])
                except ZeroDivisionError as detail:
                    pass
                    #print('Error: on this ',gene,pos, detail) 

            if (int(pos)+1)%3==2:
               
                try:
                    codon =seq[p-1:p+2]
                    ni,si=getndns(codon,codontable)
                    mutcodon= substbase(codon,1,getReverseStrand(gids[gene][pos][1]))
                    if codontable[codon]==codontable[mutcodon]:
                        ns+=(float(gids[gene][pos][4])/si)/float(gids[gene][pos][2])
                    elif codontable[codon]!=codontable[mutcodon]:
                        nd+=(float(gids[gene][pos][4])/ni)/float(gids[gene][pos][2])
                except ZeroDivisionError as detail:
                    pass
                    #print('Error: on this ',gene,pos,detail)
         
        n=np.median(np.array(n))
         
        return nd/r,ns/r,calTajD(int(n),sites,freqs)

def getReverseStrand(seq):
   fwdstrd='ACGTacgt'
   revstrd='TGCAtgca'
   trans=str.maketrans(fwdstrd,revstrd)
   return seq.translate(trans)[::-1]


def getpnpsfwd(seq,gids,codontable,gene):
   
    ns=0
    nd=0
    r= len(translateDNA(seq,codontable))
    freqs=[]
    n=[]
    sites=len(gids[gene].keys())

    if len(gids[gene].keys())>1:
    
        for pos in gids[gene].keys():
            freqs.append(float(gids[gene][pos][5]))
            n.append(float(gids[gene][pos][2]))
            
            if (int(pos)+1)%3==0:
                try:
                    codon=seq[int(pos)-2:int(pos)+1]
                    ni,si=getndns(codon,codontable)
                    mutcodon= substbase(codon,2,gids[gene][pos][1])

                    if codontable[codon]==codontable[mutcodon]:
                        ns+=(float(gids[gene][pos][4])/si)/float(gids[gene][pos][2])

                    if codontable[codon]!=codontable[mutcodon]:                    
                        nd+=(float(gids[gene][pos][4])/ni)/float(gids[gene][pos][2])

                except:
                    pass 
                    #print('Error: on this ',gene,pos) 
                
            if (int(pos)+1)%3==1:
                try:                
                    codon=seq[int(pos):int(pos)+3]
                    ni,si=getndns(codon,codontable)
                    mutcodon= substbase(codon,0,gids[gene][pos][1])
        
                    if codontable[codon]==codontable[mutcodon]:
                        ns+=(float(gids[gene][pos][4])/si)/float(gids[gene][pos][2])

                    elif codontable[codon]!=codontable[mutcodon]:
                        nd+=(float(gids[gene][pos][4])/ni)/float(gids[gene][pos][2])

                except:
                    pass
                    #print('Error: on this ',gene,pos)

            if (int(pos)+1)%3==2:
                try:
                    codon =seq[int(pos)-1:int(pos)+2]
                    ni,si=getndns(codon,codontable)
                    mutcodon= substbase(codon,1,gids[gene][pos][1])
        
                    if codontable[codon]==codontable[mutcodon]:
        
                        ns+=(float(gids[gene][pos][4])/si)/float(gids[gene][pos][2])

                    elif codontable[codon]!=codontable[mutcodon]:
                        nd+=(float(gids[gene][pos][4])/ni)/float(gids[gene][pos][2])
                except:
                    pass
                    #print('Error: on this ',gene,pos)
                    
        n=np.median(np.array(n))
        
        return nd/r,ns/r,calTajD(int(n),sites,freqs) 

def calulateNandS(seq,codontable):
  
    N=0
    S=0
    for i in range(0,len(seq)-(len(seq)%3),3):
        n,s=getndns(seq[i:i+3],codontable)
        N+=n
        S+=s
    return N,S 

def translateDNA(seq,codontable):
  
    protein=''
    for i in range(0,len(seq)-(len(seq)%3),3):
        protein+=codontable[seq[i:i+3]]
    if protein[-4:]=='STOP':
        protein=protein[:-4]
    return protein

def getgids(cdsgenes):
  
    genes={}
    with open(cdsgenes,'r') as file:
        for line in  file:
            line=line.strip().split()
            if line[0] not in genes.keys():
                genes[line[0]]={}
                genes[line[0]][line[2]]=line[3:]
            else:
                genes[line[0]][line[2]]=line[3:]
    return genes
 
def getSNPdetails(output):
  
    snps={}
    with open(output,'r') as file:
        for line in  file:
            line=line.strip().split()
            if line[0] not in snps.keys():
                snps[line[0]]=[line[1:5]]
            else:
                snps[line[0]]+=[line[1:5]]
    return snps 

def getndns(codon,codontable):   
  
    codon=codon.upper()
    bases=['A','T','G','C']
    n=0.0
    for pos,b in enumerate(codon):
        bases.remove(b)  
        for bp in bases:
            mutcod=substbase(codon,pos,bp)
            if codontable[mutcod]!=codontable[codon]:
                n+=1/3.0 
        bases=['A','T','G','C']
    return round(n,3),round(3-n,3)

def substbase(codon,pos,bp):
  
    cod=ctypes.create_unicode_buffer(codon)
    cod[pos]=bp
   
    return cod.value

def readseq(fastafile):
    sequences=[]
    seq='' 
    headers=[]

    with open(fastafile) as file:
        file1=[line.rstrip() for line in file.readlines() if len(line)>0]
        for line in file1:
            line=line.rstrip()
            if len(line)>0:
                if line[0]=='>':
                    h=line 
                    headers.append(line)
                    if seq:
                        sequences.append(seq)
                        seq=''
                else:
                    seq+=line            
    sequences.append(seq)
    return collections.OrderedDict(zip(headers,sequences))

def calTajD(n,sn,allfreqs):
   
    pi=calPi(n,allfreqs)
    theta_wat=calcThetaWats(n,sn)
    C=normalConst(n,sn)
    return (pi-theta_wat)/C

def calPi(n,allfreqs):
    c=float(n)/(n-1)
    freqsum=sum([2.0*i*(1-i) for i in allfreqs])
    return c*freqsum

def calcThetaWats(n,sn):
    """ Caclculate the Watterson's theta """
    return sn/sum([1.0/i for i in range(1,n)])


def normalConst(n,sn):
   
    a1=sum([1.0/i for i in range(1,n)])
    a2=sum([1.0/i**2 for i in range(1,n)])
    b1=float(n+1)/(3*(n-1))
    b2=float((2*(n**2+n+3)))/(9*n*(n-1))
    c1=b1-(1/a1)
    c2=b2-((n+2)/(a1*n))+a2/a1**2
    return sqrt((c1/a1)*sn + (c2/(a1**2+a2))*sn*(sn-1))

if __name__ == "__main__":
    ProcessCLI(sys.argv)
