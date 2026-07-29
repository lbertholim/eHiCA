#!/bin/bash

#cd /hihg/studies/AD/analysis/projects/jVance_functional/HiC/CaseWestern_Processed_Data/hg38/Data_analysis

request_file=$1
sample_file=$2
out_dir=$3


SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

ref="$SCRIPT_DIR/ref/hg38_DpnII_anchors_avg.bed"
chrom_size="$SCRIPT_DIR/ref/hg38.chrom.sizes"
lib="$SCRIPT_DIR/lib"

echo $SCRIPT_DIR

# for line in $(cat $SCRIPT_DIR/test.VP.loc.bed); do


tail -n +2 $request_file | while IFS=',' read -r Bait_coordinates Identifier hub_location; do

	### create directory $Bait_coordinates (chr:start-end) according to test.VP.loc.bed  
  mkdir $out_dir/$Bait_coordinates
  echo "eHiCA at" $$Bait_coordinates
  echo "1"

  #### outputs chr start end <locus> to a file .vp.bed (tab delimited)
	echo -e "$(echo $Bait_coordinates | tr ':' '\t' | tr '-' '\t')" > $out_dir/$Bait_coordinates/$Bait_coordinates.vp.bed
  echo "2"
  
  ### uses the bed file ($Bait_coordinates.vp.bed) and bedtools to extract anchors from ref file with the bed interval and outputs .anchorlist   
	bedtools intersect -wa -wb -a $ref -b $out_dir/$Bait_coordinates/$Bait_coordinates.vp.bed |  cut -f1-4 > $out_dir/$Bait_coordinates/$Bait_coordinates.anchorlist
  echo "3"

  ### Copies the anchor list into the <locus> directory for local use.
	###cp ${PWD}/$Bait_coordinates/$Bait_coordinates.anchorlist ${PWD}/$Bait_coordinates/ 

  ### Created <locus>.final - the first part of run.$Bait_coordinates.sh
  ### First Concatenates all per‑sample anchor files, creating a combined stream of anchor rows
  ### Then sorts the combined anchor entries, removes duplicate lines and outputs to union.anchorlist (to be used by heatmap.r) 
	echo cat $out_dir/$Bait_coordinates/temp"*"anchor "|" sort -u ">" $out_dir/$Bait_coordinates/union.anchorlist > $out_dir/$Bait_coordinates/$Bait_coordinates.final
  echo "4"
  ### Writes the script to run heatmap.r to <locus>.final 
	echo "Rscript $lib/heatmap.R $out_dir/$Bait_coordinates/$Bait_coordinates.anchorlist $Bait_coordinates $sample_file $out_dir"  >> $out_dir/$Bait_coordinates/$Bait_coordinates.final
  echo "5"

  ##### Runs capture.py - which reads the number of anchors in <locus>.vp.bed and chooses single‑anchor or multiple-anchor
  ###Outputs <locus>.temp(shell fragment with per‑sample capture commands) and temp.<sample>.<locus> (raw per‑sample signal)
	python $lib/capture.py \
  $out_dir/$Bait_coordinates/$Bait_coordinates.vp.bed \
  $Bait_coordinates \
  $SCRIPT_DIR \
  $sample_file \
  $out_dir;
  echo "6"

  #### Merge <locus>.temp and <locus>.final to run.$Bait_coordinates.sh file    
	cat $out_dir/$Bait_coordinates/$Bait_coordinates.temp $out_dir/$Bait_coordinates/$Bait_coordinates.final  > $out_dir/$Bait_coordinates/run.$Bait_coordinates.sh;
  echo "7"     
	echo cd $out_dir/$Bait_coordinates/ >> $out_dir/$Bait_coordinates/run.$Bait_coordinates.sh;
  echo "8"     
	echo rm temp"*" >> $out_dir/$Bait_coordinates/run.$Bait_coordinates.sh
  bash $out_dir/$Bait_coordinates/run.$Bait_coordinates.sh
  
  LOG=$out_dir/$Bait_coordinates/Hi-C_UCSC.log
  # Run VC heatmap per tissue type and export bed
  Rscript $lib/Hi-C_UCSC.r \
  $out_dir/$Bait_coordinates/Combined.captured.loop.matrix.at.$Bait_coordinates \
  $out_dir/$Bait_coordinates/$Bait_coordinates.vp.bed \
  $SCRIPT_DIR \
  $out_dir/$Bait_coordinates \
  $request_file \
  $sample_file \
  > "$LOG" 2>&1
  
  ### Run bedToBigBed
  echo "running bedTobigBed"
  bash $SCRIPT_DIR/bedTobigbed.sh \
  $out_dir/$Bait_coordinates/$Identifier/hg38 \
  $chrom_size

  echo "DONE"


done