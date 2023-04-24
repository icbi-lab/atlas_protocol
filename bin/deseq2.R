#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(conflicted)
  library("DESeq2")
  library("BiocParallel")
  library("dplyr")
  conflict_prefer("select", "dplyr", quiet=TRUE)
  library("IHW")
  library("tibble")
  library("readr")
  library("argparser", quietly = TRUE)
})

p <- arg_parser("Input Differential Expression Analysis")

p <- add_argument(p, "countData", help = "Count matrix", type = "character")
p <- add_argument(p, "colData", help = "Sample annotation matrix", type = "character")
p <- add_argument(p, "--resDir", help = "Output result directory", default = "./results")
p <- add_argument(p, "--cond_col", help = "Condition column in sample annotation", default = "group")
p <- add_argument(p, "--covariate_formula", help = "Additional covariates that will be added to the design formula", default = "")
p <- add_argument(p, "--prefix", help = "Prefix of result file", default = "test")
p <- add_argument(p, "--sample_col", help = "Sample identifier in sample annotation", default = "sample")
p <- add_argument(p, "--fdr", help = "False discovery rate", default = 0.1)
p <- add_argument(p, "--c1", help = "Contrast level 1 (perturbation)", default = "grpA")
p <- add_argument(p, "--c2", help = "Contrast level 2 (baseline)", default = "grpB")
p <- add_argument(p, "--cpus", help = "Number of cpus", default = 1)
p <- add_argument(p, "--save_workspace", help = "Save R workspace", flag = TRUE)

# Parse cmd line arguments
argv <- parse_args(p)

sampleAnno_path = argv$colData
count_mat_path = argv$countData
resDir <- argv$resDir
prefix <- argv$prefix
cond_col <- argv$cond_col
sample_col <- argv$sample_col
c1 <- argv$c1
c2 <- argv$c2
n_cpus <- argv$cpus
save_ws <- argv$save_workspace
covariate_formula <- argv$covariate_formula

# Number of cpus to be used
register(MulticoreParam(workers = n_cpus))

# Reading the Annotation sample csv file
sampleAnno <- read_csv(sampleAnno_path)
sampleAnno$sample <- sampleAnno$index
sampleAnno <- sampleAnno[,-1]
rownames(sampleAnno) <- sampleAnno$sample

# Reading the Count matrix tsv file
count_mat <- read_csv(count_mat_path)
count_mat <- count_mat |>
  select(c(gene_id, sampleAnno[[sample_col]])) |>
  column_to_rownames("gene_id") |>
  round()

# Check if names are the same
if (!all(rownames(sampleAnno) %in% colnames(count_mat))) {
  print('Row names in sample annotation and column names in count matrix are not the same')
  break
}

# Start processing
design_formula <- as.formula(paste0("~", cond_col, covariate_formula))

dds <- DESeqDataSetFromMatrix(countData = round(count_mat),
                              colData = sampleAnno,
                              design = design_formula)

## Keep only genes where  >= 10 reads per sample condition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep, ]

# run DESeq
dds.res <- DESeq(dds, parallel = (n_cpus > 1))

# Define contrast names
contrasts <- list(c(cond_col, c1, c2))
names(contrasts) <- sprintf("%s_vs_%s", c1, c2)

### IHW
# Use of IHW for p value adjustment of DESeq2 results
resIHW <- lapply(names(contrasts), function(name) {
  contrast = contrasts[[name]]
  results(dds.res, filterFun = ihw, contrast = contrast) |>
    as_tibble(rownames = "gene_id") |>
    mutate(comparison = name) |>
    arrange(pvalue)
}) |> bind_rows()


# write results to TSV files
write_tsv(resIHW, file.path(resDir, paste0(prefix, "_DESeq2_result.tsv")))

# Save R ws
if (save_ws) {
  save.image(file = file.path(resDir, paste0(prefix, "_ws.RData")))
}
