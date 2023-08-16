# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: Python [conda env:conda-2023-atlas-protocol]
#     language: python
#     name: conda-env-conda-2023-atlas-protocol-py
# ---

# %% [markdown]
# (differential_expression)=
# # Differential gene expression
#
# One of the most commonly asked questions when comparing between conditions is which genes
# are differentially expressed. In statistical terms, we ask if the mean gene expression differs between two groups more than what would be expected by random chance.
# A common fallacy here is that frequently applied statistical tests assume independence of the observations, while cells from the same sample or patient are not independent. This is also called "pseudoreplication bias" and leads to drastically inflated p-values {cite}`squairConfrontingFalseDiscoveries2021, zimmermanPracticalSolutionPseudoreplication2021`.
#
# Unfortunately, the differential testing
# functions available from the most commonly used frameworks (e.g. {func}`scanpy.tl.rank_genes_groups` in scanpy) do not account for this, allowing for qualitative comparisons at best. As an alternative, {cite:t}`squairConfrontingFalseDiscoveries2021` suggest to aggregate samples into "pseudo-bulks" and apply tools originally designed for comparing bulk RNA-seq samples.
#
# In this section, we are going to demonstrate how to generate pseudo-bulk samples by cell-type using
# [decoupler](https://github.com/saezlab/decoupler-py) {cite}`badia-i-mompelDecoupleREnsembleComputational2022` and apply DESeq2 {cite}`loveModeratedEstimationFold2014` to perform cell-type-specific comparison between conditions.
#
# :::{seealso}
#  * The [Differential gene expression](https://www.sc-best-practices.org/conditions/differential_gene_expression.html) chapter of the single-cell best practice book {cite}`heumosBestPracticesSinglecell2023`.
#  * The [decoupler pseudobulk vignette](https://decoupler-py.readthedocs.io/en/latest/notebooks/pseudobulk.html).
# :::

# %% [markdown]
# ## 1.Load libraries

# %%
# %load_ext autoreload
# %autoreload 2
import re
import subprocess as sp
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
from anndata import AnnData
from tqdm.contrib.concurrent import process_map

import atlas_protocol_scripts as aps

# %% [markdown]
# ## 2. Load input data
# 1. Load the annotated AnnData object

# %%
adata_path = "../../data/input_data_zenodo/atlas-integrated-annotated.h5ad"
adata = sc.read_h5ad(adata_path)

# %% [markdown]
# 2. Define the paths to the scripts for DESeq2. The script is shipped as part of the `atlas_protocol` repo.

# %%
deseq_script_path = "../../bin/deseq2.R"

# %% [markdown]
# 3. Create output directories

# %%
deseq_results = Path("../../data/results/differential_expression/")
# !mkdir -p {deseq_results}

# %% [markdown]
# ## 3. Generate Pseudobulk
#
# Here, we aggregate counts by biological replicate (patient) for each cell-type individually.
#
# 1. subset AnnData to the cells of interest. In this case, we only focus on cells originating from primary tumor samples and the two conditions (LUAD vs. LUSC) we are going to compare.

# %%
adata = adata[adata.obs["origin"].isin(["tumor_primary"]) & adata.obs["condition"].isin(["LUAD", "LUSC"])]
adata

# %% [markdown]
# 2. Generate pseudobulk
#
# :::{important}
# Generate pseudo-bulk based on **raw counts** and aggregate them by summing up.
# :::

# %%
# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col="sample",
    groups_col=["cell_type_coarse", "condition"],
    layer="raw_counts",  # layer where raw counts are store in adata
    mode="sum",
    min_cells=0,  # we are setting this to 0 and filter in an explicit, separate step.
    min_counts=0,
)
pdata

# %% [markdown]
# 3. Inspect the table of "pseudobulk samples"

# %%
pdata.obs

# %% [markdown]
# ## 4. Filter samples
#
# In this step, we are removing samples that consist of too few cells or have a small number of total counts in both conditions.
#
# 1. Plot relationship of counts and cells per dataset and cell-type

# %%
for col in ["dataset", "cell_type_coarse", "condition"]:
    dc.plot_psbulk_samples(pdata, groupby=col, figsize=(3, 3))

# %% [markdown]
# 2. Remove samples with less than 10 cells or less than 1000 counts

# %%
pdata = pdata[(pdata.obs["psbulk_n_cells"] >= 10) & (pdata.obs["psbulk_counts"] >= 1000)]

# %% [markdown]
# ## 5. Split pseudobulk by cell-type
#
# For the following steps, we treat each cell-type separately, therefore we split the pseudobulk object
# into a dictionary of objects per cell-type.

# %% [markdown]
# 1. Define list of cell-types

# %%
cell_types = pdata.obs["cell_type_coarse"].unique()

# %% [markdown]
# 2. Generate per-cell-type dictionary

# %%
pdata_by_cell_type = {}
for ct in cell_types:
    pb = pdata[pdata.obs["cell_type_coarse"] == ct, :].copy()
    if pb.obs["condition"].nunique() <= 1:
        print(f"Cell type {ct} does not have samples in all groups")
    else:
        pdata_by_cell_type[ct] = pb

# %% [markdown]
# ## 6. Filter genes
# Here, we filter noisy genes with a low expression value. We remove genes that do not have a certain number of total counts (`min_total_count`) and have a minimum number of counts in at least $n$ samples, where $n$ equals to the number of samples
# in the smallest group (or less, if there is a sufficient number of samples).
#
# This needs to be done for each cell-type individually since each
# cell-type may express different sets of genes.
#
# :::{seealso}
# [`filterByExpr`](https://rdrr.io/bioc/edgeR/man/filterByExpr.html) in [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html).
# :::

# %% [markdown]
# 1. Inspect, for a given cell-type, the number of genes that the thresholds (here for Epithelial cells):

# %%
dc.plot_filter_by_expr(
    pdata_by_cell_type["Epithelial cell"],
    group="condition",
    min_count=10,
    min_total_count=15,
)

# %% [markdown]
# 2. Apply filter genes for all cell-types. Here, we apply the same (default) filter to all cell-types.

# %%
for tmp_pdata in pdata_by_cell_type.values():
    dc.filter_by_expr(
        tmp_pdata,
        group="condition",
        min_count=10,
        min_total_count=15,
    )


# %% [markdown]
# ## 7. Run DESeq

# %% [markdown]
# 1. Inspect the help page of the DESeq2 script

# %%
# !{deseq_script_path} --help

# %% [markdown]
# 2. Define a function to create an output directory per cell-type


# %%
def _create_prefix(cell_type):
    ct_sanitized = re.sub("[^0-9a-zA-Z]+", "_", cell_type)
    prefix = deseq_results / "LUAD_vs_LUSC" / ct_sanitized
    prefix.mkdir(parents=True, exist_ok=True)
    return prefix


# %% [markdown]
# 2. Export pseudobulk objects to CSV files. This will generate a count matrix an a samplesheet for each cell-type.

# %%
for ct, tmp_pdata in pdata_by_cell_type.items():
    prefix = _create_prefix(ct)
    aps.io.write_deseq_tables(tmp_pdata, prefix / "samplesheet.csv", prefix / "counts.csv")


# %% [markdown]
# 3. Execute the DESeq2 script
#
# :::{important}
# Always specify at least `dataset` as a covariate when performing differential expression analysis,
# to account for batch effects.
#
# It may make sense to include other covariates (sex, tumor stage, age, ...) into the model
# depending on your dataset and research question.
# :::


# %%
def _run_deseq(cell_type: str, pseudobulk: AnnData):
    """Function: Run DESeq2 on pseudobulk of all cell types from cell_type."""
    prefix = _create_prefix(cell_type)
    # fmt: off
    deseq_cmd = [
        deseq_script_path, prefix / "counts.csv", prefix / "samplesheet.csv",
        "--cond_col", "condition",
        "--c1", "LUAD",
        "--c2", "LUSC",
        "--covariate_formula", "+ dataset + sex",
        "--resDir", prefix,
        "--prefix", prefix.stem,
    ]
    # fmt: on
    with open(prefix / "out.log", "w") as stdout:
        with open(prefix / "err.log", "w") as stderr:
            sp.run(
                deseq_cmd,
                capture_output=False,
                stdout=stdout,
                stderr=stderr,
                check=True,
            )


_ = process_map(_run_deseq, pdata_by_cell_type, pdata_by_cell_type.values())


# %% [markdown]
# 4. import results

# %%
de_results = {}
for ct in pdata_by_cell_type:
    prefix = _create_prefix(ct)
    de_results[ct] = pd.read_csv(prefix / f"{prefix.stem}_DESeq2_result.tsv", sep="\t").assign(cell_type=ct)

# %% [markdown]
# ## 8. Make volcano plots

# %% [markdown]
# 1. Convert DE results to p-value/logFC matrices as used by decoupler

# %%
logFCs, pvals = aps.tl.long_form_df_to_decoupler(pd.concat(de_results.values()), p_col="padj")

# %% [markdown]
# 2. Make a volcano plot for all cell-types of interest (e.g. B cells)

# %%
fig, ax = plt.subplots()
dc.plot_volcano(logFCs, pvals, "B cell", name="B cell", top=10, sign_thr=0.1, lFCs_thr=0.5, ax=ax)
ax.set_ylabel("-log10(FDR)")
ax.set_xlabel("log2(FC)")
fig.show()

# %% [markdown]
# 4. Normalize pseudobulk to log counts per million (logCPM)

# %% [markdown]
# ## 9. Plot genes of interest as paired dotplots
#
# Here, we are visualizing the actual expression values for each sample for all genes and cell-types of interest

# %% [markdown]
# 1. Define list of genes to visualize

# %%
genes = ["EIF1AY", "GSTM1", "XIST", "MYRF"]

# %% [markdown]
# 2. Normalize pseudobulk to counts per million (CPM)

# %%
for tmp_pdata in pdata_by_cell_type.values():
    sc.pp.normalize_total(tmp_pdata, target_sum=1e6)
    sc.pp.log1p(tmp_pdata)

# %% [markdown]
# 3. Make plot (e.g. B cells)

# %%
# Get p-values for genes of interest
tmp_pvalues = de_results["B cell"].set_index("gene_id").loc[genes, "padj"].values

aps.pl.plot_paired(
    pdata_by_cell_type["B cell"],
    groupby="condition",
    var_names=genes,
    hue="dataset",
    ylabel="log(CPM+1)",
    panel_size=(2, 4),
    pvalues=tmp_pvalues,
    pvalue_template=lambda x: f"DESeq2 FDR={x:.3f}",
)
