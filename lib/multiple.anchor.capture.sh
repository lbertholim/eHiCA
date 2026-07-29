#!/usr/bin/bash



SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

chrom_size="$PROJECT_ROOT/ref/hg38.chrom.sizes"
ref="$PROJECT_ROOT/ref/hg38_DpnII_anchors_avg.bed"
lib="$PROJECT_ROOT/lib"

datapath=$1
name=$2 #sample name
snp=$3 #bait or <locus>
anchorlist=$4 #path to anchorlist file with chr start end A_XXXX
chrID=$5
out_dir=$6

name=$name.$snp


cd $out_dir/$snp/
 
#For each anchor on the anchorlist file
# Reads all raw anchor–anchor interaction files for this chromosome and sample
for anchor in `cat $anchorlist | cut -f4 | sed 's/A_//g'`;do
	# Keeps only interactions involving the bait anchor
	cat $datapath/${chrID}*anchor |
	# Temporarily removes the A_ prefix to simplify numeric comparisons 
	grep -w A_$anchor | sed 's/A_//g' |
	## Each interaction has two anchors. Indentify non-bait anchor 
	awk -v aa=$anchor '{if($1!=aa) {print "A_"$1,$3}else{print "A_"$2,$3}}' OFS='\t' |
	## Sorts by anchor ID. Removes zero, negative and non-valid interactions. Outputs non-bait anchor and interaction strenght 
	sort -k1 | awk '{if($2>0) print $0}' >> temp.$name

done

## run normalize.r on temp.$name		
Rscript $lib/normalize.r  temp.$name

### Extract anchor IDs only
cut -f1 temp.$name > $out_dir/$snp/temp.$name.id

# Look up coordinates in reference BED, keeps only anchors present in temp.$name.id
# chr  start  end  anchor_id (temp.$name.anchor)
cat $ref | fgrep -w -f temp.$name.id | cut -f1-4  > $out_dir/$snp/temp.$name.anchor


