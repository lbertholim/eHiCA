#!/bin/env Renv

print("running Hi-C_UCSC.r")
print("Opening libraries")

library(Gviz)
library(GenomicRanges)
library (tidyr)
library(rtracklayer)
library(RColorBrewer)
options(scipen = 999)


print("Setting up paths")

args = commandArgs(trailingOnly=TRUE)
input_file=args[1]
bait_file=args[2]
script_dir=args[3]
bait_out_dir=args[4] 
request_file = args[5]
sample_file = args[6]

####Tidying combined matrix

print("Reading combined")
print(input_file)

combined<-read.table(input_file, check.names = FALSE)
combined <- cbind(coordinates = rownames(combined), combined)
rownames(combined) <- 1:nrow(combined)

combined <- combined %>% separate(coordinates, c("chr", "coordinates"), ":", remove = TRUE)
combined <- combined %>% separate(coordinates, c("start", "end"), "-", remove = TRUE)
combined$start <- as.numeric(as.character(combined$start))
combined$end <- as.numeric(as.character(combined$end))

#### reading sample file

print("Reading sample file")
sample_df<-read.table(sample_file, check.names = FALSE, header = TRUE)

#### reading bait coordinates file

print("Reading bait file")

bait<-read.table(bait_file, check.names = FALSE)

####reading and tidying request file

print("Reading request file")
request<-read.csv(request_file, check.names = FALSE, header = TRUE)

request <- request %>% separate(Bait_coordinates, c("Bait_chr", "coordinates"), ":", remove = TRUE)
request <- request %>% separate(coordinates, c("Bait_start", "Bait_end"), "-", remove = TRUE)
request$Bait_start <- as.numeric(as.character(request$Bait_start))
request$Bait_end <- as.numeric(as.character(request$Bait_end))

##### setting variables

print("Setting variables")

chr<-combined$chr[1]
from<-combined$start[1]
to<-combined$end[nrow(combined)]

bait_from<-as.numeric(bait[[2]][1])
bait_to<-as.numeric(bait[[3]][1])

request_selected<- subset(request, request$Bait_chr==chr & request$Bait_start==bait_from,)
VC_identifier<-request_selected$Identifier[1]
hub_path<-request_selected$hub_location[1]
hub_path <- gsub('"', '', hub_path)

print(paste("working on", VC_identifier, sep=" "))

print("Setting ouput path")

#### setting output file path
##S3 path for TrackDB.txt file
S3_VC_dir<-paste(hub_path, "/", VC_identifier, "/", sep="")

####Local directory for output bed files and hub txt
dir.create(file.path(paste(bait_out_dir,"/", VC_identifier,sep="")), showWarnings = FALSE)
local_track_hub_dir_id<-paste(bait_out_dir,"/", VC_identifier,"/", sep="")

dir.create(file.path(paste(local_track_hub_dir_id, "hg38", sep="")), showWarnings = FALSE)
local_track_hub_dir_id_hg38<-paste(local_track_hub_dir_id, "hg38", "/", sep="")

####Opening trackDB part2 file, hub and genomes files
trackDB_pt2<-file(paste(bait_out_dir,"/", "trackDb_pt2.txt", sep="/"), open='wt')

hub<-file(paste(local_track_hub_dir_id, "hub.txt", sep="/"), open='wt')
genomes<-file(paste(local_track_hub_dir_id, "genomes.txt", sep="/"), open='wt')


####Building gviz annotation tracks from R.environment

print("Code will use environment UCSC tracks")

load(paste0(script_dir, "/Environment.RData"))

idxTrack<- get(paste0("idxTrack_", chr))

ucscGenes<- get(paste0("UcscTrack_", chr))


axTrack <- GenomeAxisTrack()
  
aTrack <- AnnotationTrack(start = bait_from, width = bait_to-bait_from, 
                              chromosome = chr, 
                              fill="black",
                              col.line = NULL,
                              col = 0,
                              min.width=2,
                              min.height=5,
                              genome = "hg38", name = "BAIT") 



###setting colors of heatmap for each tissue/sample type
print("Building figure for each tissue/sample type")



heatmap_cat_list<-unique(sample_df$Heatmap_grouping)

for (heatmap_cat in heatmap_cat_list){
  
  print(paste("working on ", heatmap_cat, sep=""))
  
  ####Building sample heatmap
  sample_heatmap_cat<-sample_df[sample_df$Heatmap_grouping==heatmap_cat, ]
  sample_set_heatmap<-unique(sample_heatmap_cat$Sample_name)
  
  #####Building Combined heatmap
  col_names<-c("chr", "start", "end")
  col_names_heatmap<-append(col_names, sample_set_heatmap)
  combined_heatmap <- combined[, colnames(combined) %in% col_names_heatmap]

  ####Creating group list to be used to build TrackDB
  group_labels<-list()
  i<-0

  
  type_list<-unique(sample_heatmap_cat$CellType)
  
  for (type in type_list){
  
  group_list<-unique(sample_heatmap_cat$Group)
  
    for (group in group_list){
    sample_heatmap_group<-sample_heatmap_cat[sample_heatmap_cat$Group==group&sample_heatmap_cat$CellType==type, ]
    sample_set_heatmap_group<-unique(sample_heatmap_group$Sample_name)
    print(paste(type, group))
    print(sample_set_heatmap_group)
      
    group_column_values=apply(combined_heatmap[,sample_set_heatmap_group],1,mean)
    column_name<- paste(type, group, sep="_")
    combined_heatmap[column_name]<-group_column_values


    #### rbind group info do sample_df to be used by TRackDB
    sample_df_group_line<- c("NA", paste(type,group,sep="_"),group, type, heatmap_cat)
    sample_df<-rbind(sample_df, sample_df_group_line)

    ###Appending group names to a list to be used by TrackDB
    i<-i+1
    print(i)
    group_labels[[i]] <- paste0(column_name, i)
  

    
}
  
  breaks<-unique(c(seq(0,quantile(as.matrix(combined_heatmap[,-(1:3)]),.90),len=60),
                 seq(quantile(as.matrix(combined_heatmap[,-(1:3)]),.90),quantile(as.matrix(combined_heatmap[,-(1:3)]),.99),len=30),
                 seq(quantile(as.matrix(combined_heatmap[,-(1:3)]),.99),max(combined_heatmap[,-(1:3)]),len=10)))

  }
 
  

  col<-colorRampPalette(c("blue","white", "red"))(97)
  labels<-seq(1:97)
  
  combined_bin<-combined_heatmap[1:3]
  
  for (i in colnames(combined_heatmap[,-(1:3)])){
    bin_sample_column<-cut(combined_heatmap[[i]],
                           breaks = breaks,
                           include.lowest = T,
                           right = F,
                           labels=labels)
    bin_sample_column<-as.numeric(bin_sample_column)
    combined_bin[[i]]<-bin_sample_column
  }
  

 
  dTrack2 <- DataTrack(range = combined_bin[,-1], 
                       genome = "hg38", 
                       type = "heatmap", 
                       chromosome = chr, 
                       name = "heatmap",
                       ncolor = 97,
                       gradient = col)
  rownames(dTrack2@data)<-gsub("^X", "", rownames(dTrack2@data))
  

  pdf(paste(bait_out_dir, "/", VC_identifier,".",type, ".", chr,":",bait_from,"-",bait_to,".HeatMap_Color_RGB.pdf",sep =""), width=15,7)
  try(plotTracks(list(idxTrack,axTrack, ucscGenes, aTrack, dTrack2), 
                 from = from, to = to, 
                 transcriptAnnotation = "symbol",
                 showSampleNames = TRUE,
                 cex.sampleNames = 0.8,
                 col.sampleNames="black",
                 col.axis="transparent",
                 cex.title=1.0,
                 col.title="black",
                 sizes=c(1,1,3,1, 10),
                 # margin = 0,
                 # innerMargin = -1,
                 # title.width = 8.5,
                 # labelPos="below",
                 # grid=TRUE,
                 stackHeight=0.5,
                 showColorBar=F
  ))
  dev.off()

  
  library(pheatmap)
  pheatmap(t(combined_heatmap[,-(1:3)]),cluster_rows=F,cluster_cols=F,col=colorRampPalette(c("blue","white", "red"))(100),
           show_colnames=F,breaks=breaks,border="grey60",width=15,height=7,file=paste(bait_out_dir, "/", VC_identifier,".",type, ".", chr,":",bait_from,"-",bait_to,".HeatMap_Color_RGB_pheatmap.pdf",sep =""))

  
  #################### building color matrix
  
  combined_bin2<-combined_heatmap[1:3]
  
  
  col<-as.vector(colorRampPalette(c("blue","white", "red"))(97))
  col_rgb <- col2rgb(col)
  col_rgb2<- split(col_rgb, rep(1:ncol(col_rgb), each = nrow(col_rgb)))
  
  i<-0
  for (sample in colnames(combined_heatmap[,4:(ncol(combined_heatmap))])){
    print(sample)
    i<-i+1
    bin_sample_column<-cut(combined_heatmap[[sample]],
                           breaks = breaks,
                           include.lowest = T,
                           right = F,
                           labels=col_rgb2)
    
    #bin_sample_column<-as.numeric(bin_sample_column)
    ####combined all samples
    combined_bin2[[sample]]<-bin_sample_column
 
    
       
    ###per sample bed
    combined_sample<-combined_heatmap[1:3]
    #combined_sample[["name"]]<-round(combined_heatmap[[sample]]*100, 0)
    combined_sample[["name"]]<-round(combined_heatmap[[sample]], 3)
    combined_sample[["score"]]<-round(combined_heatmap[[sample]]*100, 0)
    combined_sample[["strand"]]<-"."
    combined_sample[["color_start"]]<-combined_sample$start
    combined_sample[["color_end"]]<-combined_sample$end
    bin_sample_column<-gsub("^c\\(|\\)$", "", bin_sample_column)
    bin_sample_column<-gsub("[[:space:]]", "", bin_sample_column)
    combined_sample[[sample]]<-bin_sample_column
    
    #### Replace score >100 by 100 (UCSC browser does not deal with score >100)
    combined_sample$score[combined_sample$score > 1000] <- 1000
    
    ###bait color
    condition<-combined_sample$end > bait_from & combined_sample$start < bait_to
    combined_sample[[sample]][condition]<- "0,0,0"
    
    #save.image(file="/hihg/studies/AD/analysis/projects/jVance_functional/HiC/CaseWestern_Processed_Data/hg38/scripts/lib/Environment_test.RData")

    ###retrive sample information from sample_df 
    sample_df_ind<-subset(sample_df, sample_df$Sample_name==sample)
    group<-sample_df_ind$Group[[1]]
    cellType<-sample_df_ind$CellType[[1]]
    
    #### export sample bed to file
    sample_file_path<-paste(local_track_hub_dir_id_hg38,i,"_",type,"_", chr,"_",bait_from,"-",bait_to,"-",sample,".bed",sep='')
    file=file(sample_file_path, open='wt')
    write.table(combined_sample,file,sep='\t',quote=F,row.names=F, col.names=FALSE, append=T)
    close(file)
    
    #### write trackDb.txt file
    writeLines(paste("  track ", sample, sep=''), con=trackDB_pt2)
    writeLines("  type bigBed 9", con=trackDB_pt2)
    writeLines(paste("  bigDataUrl ", i,"_",type,"_", chr,"_",bait_from,"-",bait_to,"-",sample,".bb", sep=''), con=trackDB_pt2)
    writeLines("  itemRgb On", con=trackDB_pt2)
    writeLines("  visibility dense", con=trackDB_pt2)
    if (!(sample %in% sample_set_heatmap)){
      writeLines(paste("  shortLabel ", sample, sep=''), con=trackDB_pt2)
      writeLines(paste("  longLabel ", sample, sep=''), con=trackDB_pt2)
      writeLines(paste("  subGroups type=All_samples_by_group"), con=trackDB_pt2)
      writeLines(paste("  parent", VC_identifier,"on", sep=' '), con=trackDB_pt2)
    } else {
      writeLines(paste("  shortLabel ", sample, "-",group, sep=''), con=trackDB_pt2)
      writeLines(paste("  longLabel ", sample, "-",group, sep=''), con=trackDB_pt2)
      writeLines(paste("  subGroups type=", cellType, sep=''), con=trackDB_pt2)
      writeLines(paste("  parent", VC_identifier,"off", sep=' '), con=trackDB_pt2)
    }
    writeLines("", con=trackDB_pt2)
   }
  
  write.csv(combined_bin2, file=paste(bait_out_dir, "/", VC_identifier,".",type, ".", chr,":",bait_from,"-",bait_to,".HeatMap_Color_RGB.csv",sep =""))
  write.csv(combined_heatmap, file=paste(bait_out_dir, "/", VC_identifier,".",type, ".", chr,":",bait_from,"-",bait_to,".HeatMap_Loop_strength.csv",sep =""))
}

print("Writing output files")

####write BAIT file
bait_file_path<-paste(local_track_hub_dir_id_hg38,chr,"_",bait_from,"-",bait_to,"-BAIT",".bed",sep='')
bait_row<-list(chr, bait_from, bait_to, 0, 0,".", bait_from, bait_to, "0,0,0")
file=file(bait_file_path, open='wt')
write.table(bait_row,file,sep='\t',quote=F,row.names=F, col.names=FALSE, append=T)
close(file)

#### write BAIT file to trackDb.txt file
writeLines("  track BAIT", con=trackDB_pt2)
writeLines("  type bigBed 9", con=trackDB_pt2)
writeLines(paste("  bigDataUrl ", chr,"_",bait_from,"-",bait_to,"-BAIT",".bb", sep=''), con=trackDB_pt2)
writeLines("  itemRgb On", con=trackDB_pt2)
writeLines("  visibility dense", con=trackDB_pt2)
writeLines("  shortLabel BAIT", con=trackDB_pt2)
writeLines(paste("  longLabel BAIT_", chr,"_",bait_from,"-",bait_to,sep=''), con=trackDB_pt2)
writeLines("  subGroups type=BAIT", con=trackDB_pt2)
writeLines(paste("  parent", VC_identifier,"on", sep=' '), con=trackDB_pt2)
close(trackDB_pt2)


#### write hub.txt file
writeLines(paste("hub Heatmap_eHiCA",VC_identifier,chr,bait_from,bait_to, sep="_"), con=hub)
writeLines(paste("shortLabel Heatmap eHiCA",VC_identifier, chr, bait_from,bait_to, sep=" "), con=hub)
writeLines(paste("longLabel Heatmap eHiCA - single samples with color break by type + group by tissue/cell type -", VC_identifier, sep=" "), con=hub)
writeLines("genomesFile genomes.txt", con=hub)
writeLines("email lxb2440@med.miami.edu", con=hub)
close(hub)

#### write genomes.txt file
writeLines("genome hg38", con=genomes)
writeLines("trackDb hg38/trackDb.txt", con=genomes)
close(genomes)


####Opening trackDB file and merging with trackDB part 2

celltype_list<-unique(sample_df$CellType)
celltype_list_str <- paste(paste(celltype_list, celltype_list, sep = "="), collapse = " ")

trackDB<-file(paste(local_track_hub_dir_id_hg38, "trackDb.txt", sep="/"), open='wt')
writeLines(paste("track",VC_identifier, sep=" "), con=trackDB)
writeLines("type bigBed 9", con=trackDB)
writeLines("compositeTrack on", con=trackDB)
writeLines(paste("shortLabel",VC_identifier, sep=" "), con=trackDB)
writeLines(paste("longLabel ", VC_identifier,"_", chr,"_",bait_from,"-",bait_to, sep=''), con=trackDB)
writeLines(paste("subGroup1 type Type All_samples_by_group=All_samples_by_group ", celltype_list_str, sep=''), con=trackDB)
writeLines("sortOrder type=+", con=trackDB)
writeLines("dragAndDrop on", con=trackDB)
writeLines("visibility dense", con=trackDB)
writeLines("", con=trackDB)

trackDB_pt2 <- readLines(paste(bait_out_dir,"/", "trackDb_pt2.txt", sep="/"))
writeLines(trackDB_pt2, con = trackDB)

#### write variable.txt file
main_track_hub_build_url<-paste("https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&ncbiRefSeqCurated=dense&position=",
                                chr,"%3A",sum(bait_from,-2500000),"-",sum(bait_to,2500000),"&hubUrl=",S3_VC_dir,"hub.txt",sep='')

var_data=data.frame(paste(chr, "-",bait_from,"-",bait_to, sep=""), VC_identifier, main_track_hub_build_url)
var_file=paste(bait_out_dir,"/", "track_hub_url.txt",sep='')
write.table(var_data,var_file,sep='\t',quote=F,row.names=F,col.names=F)

print(bait_out_dir)

print(paste("END OF R Hi-C_UCSC.r - ", VC_identifier, sep=""))