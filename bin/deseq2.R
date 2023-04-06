#!/usr/bin/env Rscript

'deseq2.R

Usage:
  deseq2.R <sample_sheet> <count_table> --c1=<c1> --c2=<c2> [options]
  deseq2.R --help

Arguments:
  <sample_sheet>                CSV file with the sample annotations.
  <count_table>                 TSV file with the read counts

Mandatory options:
  --c1=<c1>                     Contrast level 1 (perturbation). Needs to be contained in condition_col.
  --c2=<c2>                     Contrast level 2 (baseline). Needs to be contained in condition_col.


' -> doc


library("DESeq2")
library("docopt")
arguments = docopt

# Input and output
sampleAnnotationCSV <- arguments$sample_sheet
readCountFile <- arguments$count_table
cond_col = "group"
sample_col = "sample"
covariate_formula = ""
contrast = c(cond_col, arguments$c1, arguments$c2)

# Testdata
#sampleAnnotationCSV = "/data/projects/2023/atlas_protocol/results/differential_expression/Tumorcells_samplesheet.csv"
#readCountFile = "/data/projects/2023/atlas_protocol/results/differential_expression/Tumorcells_counts.tsv"
#result_dir ="/data/projects/2023/atlas_protocol/results/differential_expression"
#cond_col = "group"
#sample_col = "sample"
#contrast = c("group", "LUAD", "LUSC")
#covariate_formula = ""


count_mat <- as.matrix(read.csv(readCountFile,sep="\t",row.names="gene_symbol"))
count_mat <- as.data.frame(count_mat)
allSampleAnno <- read.csv(sampleAnnotationCSV, row.names=1)

colnames(count_mat) <- NULL
count_mat[,-1]= round(count_mat[,-1],0) 

design_formula <- as.formula(paste0("~", cond_col, covariate_formula))


################# Start processing
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = allSampleAnno,
                              design = design_formula)
# run DESeq
dds.res<- DESeq(dds)

### IHW

# use of IHW for p value adjustment of DESeq2 results
resIHW <- results(dds.res, filterFun=ihw, contrast=contrast)

resIHW <- as.data.frame(resIHW ) |>
  rownames_to_column(var = "gene_id") |>
  as_tibble() |>
  arrange(padj)


#### write results to TSV and XLSX files
write.csv(resIHW , "/data/projects/2023/atlas_protocol/results/differential_expression/IHWallGenes.tsv", row.names=FALSE)


