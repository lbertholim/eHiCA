args = commandArgs(trailingOnly=TRUE)
anchor=args[1]
snp=args[2]
sample_file=args[3]
out_dir=args[4]

###Reads the anchor list, extracts numeric anchor positions, keeps a BED‑style copy (vpbed) for later use
vp=read.table(anchor,sep='\t',stringsAsFactors=F)
vpbed=vp
vp=as.numeric(gsub("A_","",vp[,4]))

###sample order, sample IDs used in filenames

sample <- read.table(file.path(sample_file), sep='\t', stringsAsFactors=F, header=T)

sampleid=sample$Sample_name

###Creates a list with location of each temp.<sample>.locus (one per sample)
file=paste(out_dir, "/", snp, "/temp.", sampleid, ".", snp, sep="")

#### Ensure all anchors are included: union across samples & the bait anchor(s) - as fullid table
### Fullid is the master anchor file
id=read.table(paste(out_dir, "/", snp, "/union.anchorlist", sep=""),sep='\t',stringsAsFactors=F)
id=rbind(id,vpbed)
fullid=id[!duplicated(id[,4]),]
rownames(fullid)=fullid[,4]

###Creates data df with anchors name as first column
data=data.frame(id=unique(id[,4]))
rownames(data)=data$id

### Build Per‑sample loop matrices 
for(i in 1:length(sampleid)){

  ### reads each temp.<sample>.locus file
  df=read.table(file[i],sep='\t',stringsAsFactors=F,row.names=1)
  
  ### Finds rows in fullid whose anchor_id matches
  temppos=fullid[fullid[,4]%in%rownames(df),]
  
  ### Keeps only anchors that actually have signal in this sample
  temppos=temppos[rownames(df),]
  
  
  #### binds sample df to data
  data[,i+1]=0
  data[rownames(df),i+1]=df[,1]
}

### Defines the genomic window (range of anchors) that will be included in the combined matrix and heatmap, 
data$id=as.numeric(gsub("A_","",data$id))
start=min(vp)-400
end=max(vp)+400
if(start<0){
  start=0
}

#### Forces bait anchors to exist in the matrix

### Turns anchors ids into dataframe with same number of columns (samples) as data
### Creates a zero‑filled matrix for bait anchors
vp=data.frame(vp)
tempmat=data.frame(matrix(0,nrow(vp),ncol(data)))
colnames(tempmat)=colnames(data)

### Assign bait anchor IDs
tempmat[,1]=vp[,1]
rownames(tempmat)=paste("A_",vp[,1],sep='')


### Keeps only anchors within the ±400 anchor window around the bait
data=data[data$id%in%seq(start,end,by=1),]

### Add bait anchors if missing
data=rbind(data,tempmat)
data=data[!duplicated(data$id),]

### Order anchors by genomic position
data=data[order(data$id),]
id=data$id

### Final cleanup
### No metadata columns remain.
### Final matrix used for plotting.
data$id=NULL
colnames(data)=sampleid

plotdata=data

rownames(fullid)=fullid[,4]

### creating pos column
fullid$pos=paste(fullid[,1],paste(fullid[,2],fullid[,3],sep='-'),sep=':')
colnames(fullid)=c("chr", "start", "end", "anchor", "pos")

### Align fullid to the matrix order
fullid=fullid[rownames(plotdata),]
rownames(plotdata)=fullid$pos
outfilename=paste(out_dir, "/", snp, "/","Combined.captured.loop.matrix.at.",snp,sep='')
write.table(plotdata,outfilename,sep='\t',quote=F)