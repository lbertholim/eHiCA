#!/usr/bin/bash

echo "running single.anchor"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

chrom_size="$PROJECT_ROOT/ref/hg38.chrom.sizes"
ref="$PROJECT_ROOT/ref/hg38_DpnII_anchors_avg.bed"
lib="$PROJECT_ROOT/lib"


datapath=$1
name=$2 #sample name
snp=$3 #bait or <locus>
anchorlist=$4 #anchor
chrID=$5

echo $name

name=$name.$snp

cat $datapath/${chrID}*anchor | grep "A_$anchorlist" | head

cd $snp/

# Keeps only interactions involving the bait anchor
cat $datapath/${chrID}*anchor | 

# Temporarily removes the A_ prefix to simplify numeric comparisons 
grep -w A_$anchorlist | sed 's/A_//g' | 

## Each interaction has two anchors. Indentify non-bait anchor 
awk -v aa=$anchorlist '{if($1!=aa) {print "A_"$1,$3}else{print "A_"$2,$3}}' OFS='\t' | 

## Sorts by anchor ID. Removes zero, negative and non-valid interactions. Outputs non-bait anchor and interaction strenght
sort -k1 | awk '{if($2>0) print $0}'> temp.$name

### Extract anchor IDs
echo "4"
cut -f1 temp.$name > temp.$name.id

# Look up coordinates in reference BED, keeps only anchors present in temp.$name.id
# chr  start  end  anchor_id (temp.$name.anchor)
echo "5"
cat $ref | fgrep -w -f temp.$name.id | cut -f1-4 |sort -k4 > temp.$name.anchor

echo "end"
echo $name

