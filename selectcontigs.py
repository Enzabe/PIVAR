import string,re,sys,copy, math, collections 

def ProcessCLI(args):
    
    contigs= readseq(args[1])
    out=open(args[3], "w")
    for contig  in contigs.keys():
        if len(contigs[contig])>=int(args[2]):
            out.write(contig+"\n")
            out.write(contigs[contig]+"\n")
            #print(contig) 
            #print(contigs[contig])
            
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

if __name__ == "__main__":
    ProcessCLI(sys.argv)
