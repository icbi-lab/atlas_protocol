#!/bin/bash
# https://www.biostars.org/p/140471/#140798

usage() {
    echo "Usage: $0 -f input_file -o output_file -a annotation"
    echo "Possible options for annotation: gencode, ensembl"
    exit 1
}

while getopts ":f:o:a:" opt; do
    case ${opt} in
        f )
            file=$OPTARG
            ;;
        o )
            out=$OPTARG
            ;;
        a )
            if [ "$OPTARG" = "gencode" ]; then
                anno="gencode"
            elif [ "$OPTARG" = "ensembl" ]; then
                anno="ensembl"
            else
                echo "Invalid annotation type provided."
                usage
            fi
            ;;
        \? )
            usage
            ;;
        : )
            echo "Invalid option: $OPTARG requires an argument"
            usage
            ;;
    esac
done
shift $((OPTIND -1))

if [ -z "$file" ] || [ -z "$out" ] || [ -z "$anno" ]; then
    echo "Missing required arguments."
    usage
fi

if [ ! -f "$file" ]; then
    echo "Input file not found."
    usage
fi

if [ "${anno}" == "gencode" ]; then
    zcat "$file" | \
    awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[2]"\t"$7}' | \
    sed 's/gene_id "//' | \
    sed 's/gene_id "//' | \
    sed 's/gene_type "//'| \
    sed 's/gene_name "//' | \
    sed 's/"//g' | \
    awk 'BEGIN{FS="\t"}{split($3,a,"[:-]"); print $1"\t"$2"\t"a[1]"\t"a[2]"\t"a[3]"\t"$4"\t"$5"\t"a[3]-a[2];}' | \
    sed "1i\Geneid\tGeneSymbol\tChromosome\tStart\tEnd\tClass\tStrand\tLength" > "$out"
elif [ "${anno}" == "ensembl" ]; then
    zcat "$file" | \
    awk 'BEGIN{FS="\t"}{split($9,a,";"); if($3~"gene") print a[1]"\t"a[3]"\t"$1":"$4"-"$5"\t"a[5]"\t"$7"\t"gensub(/.*gene_version \"([^\"]+)\".*/, "\\1", "g", $9)"\t"$2}' | \
    sed 's/gene_id "//' | \
    sed 's/gene_id "//' | \
    sed 's/gene_biotype "//'| \
    sed 's/gene_name "//' | \
    sed 's/gene_biotype "//' | \
    sed 's/"//g' | sed 's/ //g' | \
    sed '1igene_id\tGeneSymbol\tChromosome\tClass\tStrand\tGeneVersion\tGeneSource' > "$out"
else
    echo "Invalid annotation type provided."
    usage
fi
