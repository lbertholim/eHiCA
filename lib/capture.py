#!/bin/python

import subprocess
import numpy as np
import pandas as pd
import sys
import os

### Helper: readin()
### Reads sample.info, skip header and assigns column names

def readin(file):
	data=pd.read_csv(file,sep="\t",names=['datapath','sample','group','CellType','Heatmap_grouping'],index_col=False,header=0)
	return (data)

###Inputs

### vpfile: bed file with bait coordinates
### vpindex: bait coordinates
### lib: lib location
vpfile=sys.argv[1]
vpindex=sys.argv[2]
script_dir=sys.argv[3]
sample_file=sys.argv[4]
out_dir=sys.argv[5]


### Files location
sampleinfo=sample_file
single=script_dir + "/lib/single.anchor.capture.sh"
mult=script_dir + "/lib/multiple.anchor.capture.sh"
anchorlist=out_dir+"/"+vpindex+"/"+vpindex+".anchorlist"

### reading <locus>.anchorlist bed file(chr	start	end	A_XXXX), which are anchors overlapping the locus

size=pd.read_csv(anchorlist,sep='\t',index_col=False,header=0)


#### <locus>.temp table contruction
## each row = one sample, and columns correspond to future shell arguments
sample=readin(sampleinfo)
sample["vpfile"]=vpfile
sample["vpindex"]=vpindex
sample["out_dir"]=out_dir
#order=["vpfile",'datapath','sample',"vpindex"]
#sample=sample[order]
if len(size)==0:
    anchors=list(size)[3].replace("A_","")
    sample["capture"]=single
    sample['anchors']=anchors
    sample["chr"]=list(size)[0]
else:
    sample["capture"]=mult
    sample['anchors']=anchorlist
    sample["chr"]=size.iloc[0,0]

sample['com']="/usr/bin/bash"
order=["com","capture",'datapath','sample',"vpindex","anchors","chr", "out_dir"]
sample=sample[order]

sample.to_csv(out_dir+"/"+vpindex+"/"+vpindex+".temp",sep='\t',header=False,index=False)


