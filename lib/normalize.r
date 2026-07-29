#!/bin/env Rscript

#Collapses multiple rows per anchor into a single value by taking the mean interaction strength.

args = commandArgs(trailingOnly=TRUE)

# takes temp.<sample>.<locus>
file=args[1]

# Reads the 2‑column tab‑separated file
df=read.table(file,sep='\t',stringsAsFactors=F)

data=aggregate(df[,2],by=list(df[,1]),mean)

write.table(data,file,sep='\t',quote=F,row.names=F,col.names=F)

