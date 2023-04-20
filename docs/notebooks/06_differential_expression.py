# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:conda-2023-atlas-protocol]
#     language: python
#     name: conda-env-conda-2023-atlas-protocol-py
# ---

# %% [markdown]
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
import glob
import os
import subprocess as sp

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from IPython.display import display

from atlas_protocol_scripts.pl import plot_paired

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
deseq_results = "../../data/results/differential_expression/"
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
pdata_by_cell_type = {
    cell_type: pdata[pdata.obs["cell_type_coarse"] == cell_type, :].copy() for cell_type in cell_types
}

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
for cell_type in cell_types:
    dc.filter_by_expr(
        pdata_by_cell_type[cell_type],
        group="condition",
        min_count=10,
        min_total_count=15,
    )

# %% [markdown]
# 4. Normalize pseudobulk to log counts per million (logCPM)

# %%
sc.pp.normalize_total(pdata, target_sum=1e6)
sc.pp.log1p(pdata)

# %% [markdown]
# ## 6. Create dictionary of adatas subsetted by cell type

# %%
adata_dict = {}
for name in cell_type:
    name_ad = name.replace(" ", "_")
    adata_name = f"{name_ad}_adata"
    adata_dict[adata_name] = adata[adata.obs["cell_type_coarse"].isin([name])]


# %% [markdown]
# ##  Function: Run DESeq2 on pseudobulk of all cell types from cell_type


# %%
def run_deseq(count_table, sample_sheet, deseq_prefix, contrast, deseq_resdir):
    """Function: Run DESeq2 on pseudobulk of all cell types from cell_type."""
    os.makedirs(deseq_resdir, exist_ok=True)

    # fmt: off
    deseq_cmd = [
        deseq_script_path, count_table, sample_sheet,
        "--cond_col", "condition",
        "--c1", contrast[0],
        "--c2", contrast[1],
        "--resDir", deseq_resdir,
        "--prefix", deseq_prefix,
        "--cpus", 1,
    ]
    # fmt: on
    with open(deseq_resdir + "/" + deseq_prefix + ".log", "w") as stdout:
        with open(deseq_resdir + "/" + deseq_prefix + ".err", "w") as stderr:
            sp.run(
                deseq_cmd,
                capture_output=False,
                stdout=stdout,
                stderr=stderr,
                check=True,
            )


# %% [markdown]
# ## Function: Homogenize all sample ids


# %%
def fix_sample_ids(pb):
    """Homogenize all sample ids."""
    repl = {}
    for k, v in dict(zip(pb.obs["condition"].index, "_" + pb.obs["condition"].values)).items():
        repl[k] = k.replace(v, "")

    return repl


# %% [markdown]
# ## Funciton: Save results from pseudobulk (samplesheet and counts) for all cell types from cell_type


# %%
def save_pseudobulk(pb, samplesheet_filename, counts_filename):
    """Save results from pseudobulk (samplesheet and counts) for all cell types from cell_type."""
    samplesheet = pb.obs.copy()
    samplesheet.reset_index(inplace=True)
    sample_ids_repl = fix_sample_ids(pb)
    bulk_df = pb.to_df().T.rename(columns=sample_ids_repl)
    bulk_df = pb.to_df().T
    bulk_df.index.name = "gene_id"
    samplesheet.to_csv(samplesheet_filename, index=False)
    bulk_df.to_csv(counts_filename)


# %% [markdown]
# ## 7. Create pseudobulk for each celltype using the coarse cell type annotation

# %%
for ct, tmp_ad in adata_dict.items():
    pb = dc.get_pseudobulk(
        tmp_ad,
        sample_col="sample",
        groups_col="condition",
        layer="raw_counts",
        mode="sum",
        min_prop=0.05,
        min_cells=10,
        min_counts=1000,
        min_smpls=2,
    )
    if pb.obs["condition"].nunique() <= 1:
        print(f"Cell type {ct} does not have enough replicates per group")
    else:
        contrast = ["LUSC", "LUAD"]
        contrast_str = f"{contrast[0]}_vs_{contrast[1]}"
        deseq_resdir = f"{deseq_results}/{contrast_str}"

        ct = ct.replace(" ", "_")
        ct_fname = ct.replace("/", "_")
        deseq_prefix = f"{contrast_str}_{ct_fname}"

        sample_sheet = f"{deseq_results}/{deseq_prefix}.samplesheet.csv"
        count_table = f"{deseq_results}/{deseq_prefix}.counts.csv"

        save_pseudobulk(pb, sample_sheet, count_table)
        run_deseq(count_table, sample_sheet, deseq_prefix, contrast, deseq_resdir)


# %%
# Contrasts
contrasts = [
    {"name": "LUSC_vs_LUAD", "condition": "LUSC", "reference": "LUAD"},
]

# %%
# Cell type name without space for file name
cell_type_fn = [j.replace(" ", "_").replace("/", "_") for i, j in enumerate(cell_type)]

# %%
# Get results from folder
for contrast in contrasts:
    de_res = {}

    for ct in cell_type_fn:
        csv_file = glob.glob(deseq_resdir + "/" + contrast["name"] + "_" + ct + "_adata_DESeq2_result.tsv")
        if len(csv_file) > 0:
            res_df = pd.read_csv(csv_file[0], sep="\t")
            res_df = res_df.set_index(["gene_id"])
            ct = ct.replace("_", " ")
            # Register cell type results
            de_res[ct] = res_df
    contrast["de_res"] = de_res

# %%
# Concat and build the log2FoldChange change matrix
for contrast in contrasts:
    lfc_mat = (
        pd.concat(
            [
                res.loc[:, ["log2FoldChange"]].rename(columns={"log2FoldChange": ct})
                for ct, res in contrast["de_res"].items()
            ],
            axis=1,
            sort=True,
        )
        .fillna(0)
        .T
    )
    contrast["lfc_mat"] = lfc_mat
    display(lfc_mat)

# %%
# Concat and build the fdr matrix
for contrast in contrasts:
    fdr_mat = (
        pd.concat(
            [res.loc[:, ["padj"]].rename(columns={"padj": ct}) for ct, res in contrast["de_res"].items()],
            axis=1,
            sort=True,
        )
        .fillna(1)
        .T
    )
    contrast["fdr_mat"] = fdr_mat
    display(fdr_mat)

# %%
# Concat and build the stat matrix
for contrast in contrasts:
    stat_mat = (
        pd.concat(
            [res.loc[:, ["stat"]].rename(columns={"stat": ct}) for ct, res in contrast["de_res"].items()],
            axis=1,
            sort=True,
        )
        .fillna(0)
        .T
    )
    contrast["stat_mat"] = stat_mat
    display(stat_mat)

# %% [markdown]
# ##  Volcano plots of expression
#
# Genereate volcano plots for each of the cell types on the DE genes with significant activity differences using the DESeq2 log2foldChange and padj values.

# %%
for contrast in contrasts:
    # Extract logFCs and pvals
    logFCs = contrast["lfc_mat"]
    pvals = contrast["fdr_mat"]

    n_sig = 0
    sig_tf = {}
    cell_type = list(logFCs.index)
    ct_dict = {}
    ct_dict["cell_type"] = []
    for ct in cell_type:
        ct_dict["cell_type"].append({"ct": ct})

    # generate a volcano plot panel for each cell type
    for ct in ct_dict.keys():
        n_sig = len(ct_dict[ct])

        # Calculate nrows based on ncol
        ncols = 4 if n_sig >= 4 else n_sig
        nrows = int(np.ceil(n_sig / ncols))
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 4, nrows * 4))
        empty_axs = iter(axs.flatten())

        for ct, ax in zip(cell_type, empty_axs):
            dc.plot_volcano(
                logFCs,
                pvals,
                ct,
                name=ct,
                top=10,
                sign_thr=0.1,
                lFCs_thr=0.5,
                return_fig=False,
                ax=ax,
            )
        for ax in empty_axs:
            ax.set_visible(False)

        plt.tight_layout()
        plt.show()

# %% [markdown]
# ##  Pairwise expression plot.
#
# Visualize the top n (e.g 5) genes for each cell type, colored by dataset.

# %%
for ct in cell_type:
    print("\t\t\t" + ct)
    plot_paired(
        pdata,
        "condition",
        var_names=list(contrast["de_res"][ct].index)[0:5],
        n_cols=5,
        panel_size=(2, 4),
        show_legend=True,
        hue="dataset",
        size=10,
        ylabel="normalized counts",
    )

# %%
