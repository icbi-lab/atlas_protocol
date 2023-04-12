#!/bin/bash
# https://www.biostars.org/p/140471/#140798
zcat $@ | \
awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' | \
sed 's/gene_id "//' | \
sed 's/gene_id "//' | \
sed 's/gene_type "//'| \
sed 's/gene_name "//' | \
sed 's/"//g' | \
awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | \
sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength"
