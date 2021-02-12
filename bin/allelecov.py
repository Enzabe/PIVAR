import string,re,sys,copy, math, collections, copy 

def ProcessCLI(args):
    readfile(args[1],readhsout(args[2]),args[3])
    #out=open(args[3],"w")


def readfile(allelefile,alleles,outfile):
    out=open(outfile,"w")
    with open(allelefile,'r') as file:

        for line in file:
            bases={}
            words=line.strip().split(None)
            contig = words[0]
            pos = words[1]
            ref = words[2]
            mut=words[3].split("/")
            depth = words[4]
            Aa = float(words[5])
            a = float(words[6])
            Ca = float(words[7])
            c = float(words[8])
            Ga = float(words[9])
            g = float(words[10])
            Ta = float(words[11])
            t = float(words[12])
            bases['A']=Aa+a 
            bases['C']=Ca+c 
            bases['G']=Ga+g 
            bases['T']=Ta+t
        
            if len(mut)==1 and bases[ref]>0:
                try:
                    out.write('\t'.join ((contig,pos,pos,ref,mut[0],depth,str(bases[ref]),str(bases[mut[0]]),str( alleles[(contig,pos)])))+"\n")
                    #print('\t'.join ((contig,pos,pos,ref,mut[0],depth,str(bases[ref]),str(bases[mut[0]]),str( alleles[(contig,pos)]))))
                except :
                    pass 
            if len(mut)>1:
                if ref in mut:
                    mut2= copy.copy(mut)
                    mut2.remove(ref)
                    if len(mut2)==1 and bases[ref]>0:
                        try:
                            out.write('\t'.join ((contig,pos,pos,ref,mut2[0],depth,str(bases[ref]),str(bases[mut2[0]]),str( alleles[(contig,pos)])))+"\n")
                            #print('\t'.join ((contig,pos,pos,ref,mut2[0],depth,str(bases[ref]),str(bases[mut2[0]]),str( alleles[(contig,pos)]))))
                        except:
                            pass 

def readhsout(hsoutfile):
    
    with open(hsoutfile,'r') as file:
        alleles={}
        for line in file :
            line=line.strip().split(None)
            alleles[(line[0],line[1])]=line[2]
    return alleles              

if __name__ == "__main__":
    ProcessCLI(sys.argv)
