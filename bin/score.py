#!/usr/bin/env python
# This will load in a BAMCov output file and look at the coverage and mean depth
import sys, re, os
import numpy as np

if len(sys.argv) != 2: 
    print ("I need a file name")    
    exit()
results_file = sys.argv[1]

def getCoverageAndMeanDepth(filename):
    lines = [re.split('\s+',thisline.rstrip()) for thisline in open(filename)]
    covdep = np.array([[x[5],x[6]] for x in lines], dtype=np.double)
    if covdep.shape[0] < 2:
        return np.concatenate(([0.], covdep[0], covdep[0], covdep[0])) 
    else:
        ref  = covdep[-1,:]   # last row is true reference values
        infr = covdep[:-1,:]  # other rows are inferred by refseq_masher
        return np.concatenate(([infr.shape[0]], np.min(ref-infr,0) , np.mean(ref-infr,0), np.max(ref-infr,0))) # axis 0: over columns


differences = getCoverageAndMeanDepth(results_file)
for dif in differences:
    print (dif, end = "\t")

print () # newline
