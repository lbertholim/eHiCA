#!/bin/bash

bed_path=$1
chrom_sizes=$2

cd $bed_path

for file in *.bed; do
    [ -e "$file" ] || { echo "No .bed files found in $bed_path"; break; }

    sorted_file="${file%.bed}_sorted.bed"

    sort -k1,1 -k2,2n "$file" > "$sorted_file"
    bedToBigBed "$sorted_file" "$chrom_sizes" "${file%.bed}.bb"

    rm "$sorted_file"
    rm "$file"
done




