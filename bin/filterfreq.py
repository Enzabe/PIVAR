import string,re,sys,copy,collections
import ctypes
from math import * 
import numpy as np  

def ProcessCLI(args):

    filterfreq(args[1],args[2])
    #out=open(args[2],"w")

def filterfreq(output,outfile):
    out=open(outfile,"w")
    with open(output,'r') as file:
        for cfreq in  file:
            line=cfreq.strip().split()
            if '1.0000'  not in line[5:]:
                out.write(cfreq.strip()+"\n") 
                
if __name__ == "__main__":
    ProcessCLI(sys.argv)
