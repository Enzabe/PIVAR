from __future__ import division
import string,re,sys,copy, math
import collections, itertools

def ProcessCLI(args):

    processMPIleup(args[1],args[2])

def processMPIleup(mpileupfile,allelefile):

    w=open(allelefile+".alleles","w")
    freqs = open(allelefile+".freqs","w")

    with open(mpileupfile,"r") as f:

        for line in f:
            c=line
            line=line.strip().split()
            pileup,INDEL,bases=_processpileup(line)

            if int(line[3])>9 and len(set(bases.upper()))>=2 and len(INDEL)<1:
                pbases=countbases(bases)

                if pbases['*']==0:
                    dbases,calls,ref,raf,nonRefAllele, nraf= _getAlleleFreq(line,pbases)

                    if nonRefAllele is not None  and dbases[nonRefAllele]>2:
                        maf=min(raf,nraf)
                        w.write(line[0]+"\t"+ line[1]+"\t"+line[2]+"\t"+nonRefAllele+"\t"+line[3]+"\t"+str(pbases['A'])+"\t"+str(pbases['a'])+"\t" + str(pbases['C'])+"\t"+str(pbases['c'])+"\t"+str(pbases['G'])+"\t"+str(pbases['g'])+"\t"+str(pbases['T'])+"\t"+str(pbases['t'])+"\n")
                        freqs.write(line[0] + "\t" + line[1] + "\t" + ref + "\t" + nonRefAllele + "\t" + str(dbases[ref]+dbases[nonRefAllele]) + "\t" + "{0:.4f}".format(raf) + "\t" + "{0:.4f}".format(nraf) + "\t"\
 + "{0:.4f}".format(maf) + "\n")

def findSNPs(line,pbases):
    
    dbases,calls,ref,raf,nonRefAllele, nraf= _getAlleleFreq(line,pbases)

    if nonRefAllele is not None  and dbases[nonRefAllele]>2:
        return line[:4],calls,ref,raf,nonRefAllele,nraf,dbases[ref],dbases[nonRefAllele]
    else: pass 

def countbases(seq):

    bases={'A':0.0,'a':0.0,'C':0.0,'c':0.0,'G':0.0,'g':0.0,'T':0.0,'t':0.0,'*':0.0}
    for c in seq:
        try:
            bases[c]+=1

        except KeyError as error:
            pass

    return bases


def _procLine(line,pileup):

    reffwd=line[2].upper()
    refrev=reffwd.lower()
    bases=pileup.replace('.',reffwd)
    bases=bases.replace(',',refrev)
    bases=bases.replace('$','')
    pbases=countbases(bases)
    return pbases 

def _processpileup(line):

    reffwd=line[2].upper()
    refrev=reffwd.lower()

    newpileup=re.sub(r"\^[\x00-\x7F]","",line[4])
    newmpileup=re.sub(r"[\+|-][\d+][\w+]","",newpileup).replace('$','')
    INDEL=re.findall(r"[\+|-][\d+][\w+]",line[4])
    
    bases=newpileup.replace('.',reffwd)
    bases=bases.replace(',',refrev)
    bases=bases.replace('$','')

    return newmpileup,INDEL,bases 

  
def _getAlleleFreq(line,pbases):

    ref=line[2].upper()
    contig = line[0]
    pos = line[1]
    depth = float(line[3])

    dbases = {}

    dbases['A'] = pbases['A']+pbases['a']+0.0
    dbases['C'] = pbases['C']+pbases['c']+0.0
    dbases['G'] = pbases['G']+pbases['g']+0.0
    dbases['T'] = pbases['T']+pbases['t']+0.0
    
    ab = 0 if (dbases['A']) == 0 else min((pbases['A']+0.0)/(dbases['A']), (pbases['a']+0.0)/(dbases['A']))
    cb = 0 if (dbases['C']) == 0 else min((pbases['C']+0.0)/(dbases['C']), (pbases['c']+0.0)/(dbases['C']))
    gb = 0 if (dbases['G']) == 0 else min((pbases['G']+0.0)/(dbases['G']), (pbases['g']+0.0)/(dbases['G']))
    tb = 0 if (dbases['T']) == 0 else min((pbases['T']+0.0)/(dbases['T']), (pbases['t']+0.0)/(dbases['T']))

    calls = []
    if ab > 0.25:
        calls.append('A')
    if cb > 0.25:
        calls.append('C')
    if gb > 0.25:
        calls.append('G')
    if tb > 0.25:
        calls.append('T')
        
    nonRefAllele, nraf=_findnonrefallele(ref,dbases,calls,pbases)
    raf=0.0
    if len(set(calls))>1 and ref in calls:
        try:
        
            raf=(dbases[ref]+0.0)/(dbases[nonRefAllele]+dbases[ref])
        except KeyError:
            print(calls,ref,nonRefAllele,ab,cb,gb,tb)


    return dbases,calls,ref,raf,nonRefAllele, nraf 

def _findnonrefallele(ref,dbases,calls,pbases):

    nrefCalls = calls[:]
    maxvaf = 0.0
    maxallele =''
    nraf =0.0
    nonRefAllele=''
    
    if len(calls)>1 and ref in nrefCalls:
        
        nrefCalls.remove(ref)
        for nonRefAllele in nrefCalls:
            
            vaf = (dbases[nonRefAllele]+0.0)/(dbases[nonRefAllele] + dbases[ref])
        
            if vaf > maxvaf:
                
                maxvaf = vaf
                maxallele = nonRefAllele
                nraf = maxvaf
                nonRefAllele = maxallele

    if len(calls)<2 or ref not in calls:
        
        nonRefAllele=None
        nraf=None 

    return nonRefAllele, nraf 

if __name__ == "__main__":
    ProcessCLI(sys.argv)



