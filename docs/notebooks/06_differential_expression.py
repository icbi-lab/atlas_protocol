# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:.conda-atlas_protocol_23]
#     language: python
#     name: conda-env-.conda-atlas_protocol_23-py
# ---

# %% [markdown]
# # Differential gene expression analysis per cell-type between conditions
#
# This notebookios is focused in scRNA-seq gene expression readouts and their differences in magnitude and significance between the condition of interest and a reference.
#
# In this notebook we use:
#
# - decoupler {cite}`Badia-i-Mompel2022` a tool that contains different statistical methods to extract biological activities from omics data using prior knowledge.
# - deseq2  {cite}`Love MI, Huber W, Anders S (2014).` a tools that uses the negative binomial distribution to model count data from high-throughput sequencing assays and evaluate the relationship between variance and mean, while also testing for differential expression.
#

# %% [markdown]
# ## 1.Load libraries

# %%
import glob
import os
import subprocess as sp

# Libraries for visualization
import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from IPython.display import display

from atlas_protocol_scripts.pl import plot_paired

# set PATH env variable to conda env for specific R version.
# To use [DESeq2, R version "4.2" required](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
path_to_R = "/usr/local/bioinf/R/R-4.2.3/bin/"
os.environ["PATH"] = path_to_R + os.pathsep + os.environ["PATH"]

cpus = 6

# %% [markdown]
# ## 2. Load input data
# * adata_path: Path to anndata file
# * deseq: Path to deseq2 script
# * deseq_results: Path to results directory.
#

# %%
adata_path = "/data/projects/2023/atlas_protocol/input_data_zenodo/atlas-integrated-annotated.h5ad"
deseq = "../../bin/deseq2.R"
deseq_results = "/data/projects/2023/atlas_protocol/results/differential_expression/deseq_resdir"

# %%
adata = sc.read_h5ad(adata_path)

# Subset adata
adata = adata[adata.obs["origin"].isin(["tumor_primary"])]
adata = adata[adata.obs["condition"].isin(["LUAD", "LUSC"])]

# %%
adata

# %% [markdown]
# ## 3. Pseudobulk
# Get pseudobulk for entire adata

# %%
# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col="sample",
    groups_col="condition",
    layer="raw_counts",  # layer where raw counts are store in adata
    mode="sum",
    min_cells=0,
    min_counts=0,
)
pdata

# %% [markdown]
# Normalize pseudobulk

# %%
pdata.layers["counts"] = pdata.X.copy()
sc.pp.normalize_total(pdata, target_sum=1e6)

# %%
pdata

# %% [markdown]
# ## 4. Quality control plot
#
# From generated profile for each dataset

# %%
dc.plot_psbulk_samples(pdata, groupby=["dataset", "condition"], figsize=(13, 5))

# %% [markdown]
# Convention to filter low quality samples:
#
#     - Number of cells --> min_cells (genes  minimum total number of reads across sample)
#     - Number of counts --> min_counts (genes minimum number of counts in a number of samples)
#
# Check the frequency of genes (features) vs n. of samples and total sum of counts with  ```plot_filter_by_expr```
#

# %%
dc.plot_filter_by_expr(pdata, group="condition", min_count=10, min_total_count=15)

# %% [markdown]
# ##  5. Define cell type to use
#
# Specifiy for which cell type annotation level we want to run the differential expression analyses

# %%
cell_type = list(adata.obs["cell_type_coarse"].unique())

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

    deseq_cmd = [
        deseq,
        count_table,
        sample_sheet,
        "--cond_col",
        "condition",
        "--c1",
        contrast[0],
        "--c2",
        contrast[1],
        "--resDir",
        deseq_resdir,
        "--prefix",
        deseq_prefix,
        "--cpus",
        str(cpus),
        "--save_workspace",
    ]
    with open(deseq_resdir + "/" + deseq_prefix + ".log", "w") as stdout:
        with open(deseq_resdir + "/" + deseq_prefix + ".err", "w") as stderr:
            sp.run(deseq_cmd, capture_output=False, stdout=stdout, stderr=stderr, check=True)


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
