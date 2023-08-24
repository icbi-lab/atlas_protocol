#!/bin/bash
# Modified from https://www.biostars.org/p/140471/#140798

# Check for correct number of arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input_file> <output_file> <gtf_source>"
  exit 1
fi

# Input and output file names
file="$1"
out="$2"
gtf_source="$3"

# Perform the zcat, awk, and sed commands based on annotation gtf source, and save the output to the specified output file
if [ "$gtf_source" == "gencode" ]; then
  zcat "$file" | \
  awk 'BEGIN{FS="\t"; OFS=","}{split($9,a,";"); if($3~"gene") print a[1], gensub(/^ +| +$/,"", "g", a[3]), $1, $4, $5, gensub(/^ +| +$/,"", "g", a[2]), $7, $5-$4;}' | \
  sed 's/gene_id "//' | \
  sed 's/gene_id "//' | \
  sed 's/gene_type "//'| \
  sed 's/gene_name "//' | \
  sed 's/"//g' | \
  sed "1i\Geneid,GeneSymbol,Chromosome,Start,End,Class,Strand,Length" > "$out"
elif [ "$gtf_source" == "ensembl" ]; then
  zcat "$file" | \
  grep -v '^#' | \
  awk '$3 == "gene"' | \
  cut -f 1,4-5,7,9 | \
  awk -F'\t' -v OFS='\t' '{ split($NF, a, /[;\"]+/); $NF=""; for (i=2; i<=length(a); i+=2) $(NF+1) = a[i]; if ($6 == "havana" || $6 == "ensembl") $6 = ""; $6 = $6 "." $7; $7 = ""; } 1' | \
  awk -F'\t' -v OFS='\t' '{ print }' | \
  cut -f 1-4,6,8,9,10 | \
  sed 's/\t/,/g' | \
  awk -F',' -v OFS=',' 'BEGIN{print "gene_id","gene_name","chromosome","start","end","gene_biotype","gene_source","strand","length"} { gsub("havana", "", $6); gsub("ensembl", "", $6); print $5, $6, $1, $2, $3, $8, $7, $4, $3-$2+1}' > "$out"
else
  echo "Invalid gtf source. Please specify 'gencode' or 'ensembl'."
  exit 1
fi
