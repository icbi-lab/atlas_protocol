#!/usr/bin/env Rscript

library("DESeq2")
#library("docopt")
#arguments = docopt(doc, version = "0.1")

# Testdata
# Example1
sampleAnnotationCSV = "/data/projects/2023/atlas_protocol/results/differential_expression/Tumorcells_samplesheet.csv"
readCountFile = "/data/projects/2023/atlas_protocol/results/differential_expression/Tumorcells_counts.tsv"
result_dir ="/data/projects/2023/atlas_protocol/results/differential_expression"
cond_col = "group"
sample_col = "sample"
contrast = c("group", "LUAD", "LUSC")
covariate_formula = ""



# Input and output
#sampleAnnotationCSV <- arguments$sample_sheet
#readCountFile <- arguments$count_table
#results_dir = arguments$result_dir
#contrast = c(cond_col, arguments$c1, arguments$c2)
#covariate_formula = arguments$covariate_formula

count_mat <- read_tsv(readCountFile)
allSampleAnno <- read_csv(sampleAnnotationCSV)

count_mat[,-1]= round(count_mat[,-1],0)

#sampleAnno <- allSampleAnno %>%
#  filter(get(cond_col) %in% contrast[2:3])

design_formula <- as.formula(paste0("~", cond_col, covariate_formula))


################# Start processing
dds <- DESeqDataSetFromMatrix(countData = count_mat[,-1],
                              colData = allSampleAnno,
                              design = design_formula)
# Set the reference to the contrast level 2 (baseline) given by the --c2 option
dds[[cond_col]] = relevel( dds[[cond_col]], contrast[[3]])


# run DESeq
dds <- DESeq(dds)

### IHW

# use of IHW for p value adjustment of DESeq2 results
resIHW <- results(dds, filterFun=ihw, contrast=contrast)

resIHW <- as.data.frame(resIHW ) |>
  rownames_to_column(var = "gene_id") |>
  as_tibble() |>
  arrange(padj)

#### write results to TSV and XLSX files
write.csv(resIHW , "/data/projects/2023/atlas_protocol/results/differential_expression/IHWallGenes.tsv", row.names=TRUE)
