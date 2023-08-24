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
# # Quality control and filtering
#
# For a meaningful downstream analysis, it is necessary to distinguish empty and low quality
# droplets (or wells) from *bona fide* cells. This can be achieved by investigating the number of detected
# UMIs and genes per cell, as well as the fraction of mitochondrial reads. A low number of genes
# and counts per cell can indicate empty droplets or microwells, respectively. A high fraction of
# mitochondrial reads, on the other hand, may indicate ruptured cells that lost most of their cytoplasmic RNA
# having retained only their mitochondria {cite}`lueckenCurrentBestPractices2019,lunStepbystepWorkflowLowlevel2016`.
# The metrics need to be considered jointly, as a high mitochondrial content *per se* could also be indicative of respiratory processes
# being upregulated in the cell, conveying a meaningful biological signal.
#
# Appropriate cutoffs are commonly determined by plotting the distributions of the quality metrics across all cells in a
# sample or dataset, and visually determining breakpoints between “signal” and “noise” distributions
#  {cite}`lueckenCurrentBestPractices2019`. This is the strategy we applied for building the lung cancer atlas in {cite}`salcherHighresolutionSinglecellAtlas2022a`.
# However, since these thresholds need to be determined for each dataset -- or, ideally, sample -- independently, this process requires a lot of manual interventions.
# Therefore, here we demonstrate automated removal of outliers based on median absolute deviation as suggested in {cite:t}`germainPipeCompGeneralFramework2020` and {cite:t}`heumosBestPracticesSinglecell2023`.
#
# :::{important}
# **datasets vs. studies**
#
# By “study”, we refer to a scientific publication, while with “dataset”,
# we refer to a set of samples that was generated using the same sequencing platform and processed in the same way.
# One study may contain one or multiple datasets. Datasets must be processed independently.
# :::
#
# :::{seealso}
# The [Quality Control](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html) chapter of the single-cell best practice book {cite}`heumosBestPracticesSinglecell2023`.
# :::
#

# %% [markdown]
# ## 1. Load the required libaries

# %%
import anndata
import numpy as np
import scanpy as sc

import atlas_protocol_scripts as aps

# %% [markdown]
# ## 2. Load input data
#
# TODO: use merged object instead of quick-and-dirty merge

# %%
out_dir = "../../data/results/qc/"
# !mkdir -p {out_dir}

# %%
DATASETS = {
    "maynard_2020": "../../data/input_data_raw/maynard2020.h5ad",
    "lambrechts_2018": "../../data/input_data_raw/lambrechts_2018_luad_6653.h5ad",
    "ukim-v": "../../data/input_data_raw/ukim_v_batch1.h5ad",
}
datasets = {dataset_id: sc.read_h5ad(path) for dataset_id, path in DATASETS.items()}

# %%
for ad in datasets.values():
    ad.var_names_make_unique()

# %%
adata = anndata.concat(datasets, join="inner", index_unique="_")

# %%
adata.obs["dataset"] = adata.obs_names.str.extract(r"_(.*)$", expand=False).str.replace("^\\d+_", "", regex=True)
adata.obs["patient"] = adata.obs["patient"].astype(str)

# %% [markdown]
# :::{note}
# **Ambient RNA removal**
#
# Both droplet and microwell based sequencing are subject to ambient RNA contamination. These are
# RNA molecules that are uniformly present in the cell suspension and may, for instance, originate
# from dead cells. As a consequence, ambient RNA molecules are profiled together with the
# cell-specific RNA in each droplet or well. This introduces a bias in the data that may hamper
# downstream analysis. For instance, if droplets are contaminated with RNA of cell-lineage markers,
# their expression may show up in cell-types that are known not to express these genes. Ambient
# RNA can be removed computationally by tools such as SoupX {cite}`youngSoupXRemovesAmbient2020`,
# DecontX {cite}`yangDecontaminationAmbientRNA2020`, CellBender {cite}`flemingUnsupervisedRemovalSystematic2022`,
# and SCAR {cite}`shengProbabilisticMachineLearning2022` that statistically model the measured UMI counts in each cell as the mixture
# of cell-specific and cell-free RNA.
#
# Ambient RNA removal methods typically require unfiltered UMI counts as input, which is not routinely
# available for publicly avilable dataset. For this reason, we did not perform ambient RNA removal {cite}`salcherHighresolutionSinglecellAtlas2022a`,
# and are not showing it as part of this tutorial.
# :::

# %% [markdown]
# ## 3. Calculate QC metrics
#
# 1. Label genes by common QC categories, such as mitochondrial genes, ribosomal genes and hemoglobin genes

# %%
# mitochondrial genes
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes.
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

# %% [markdown]
# 2. Add per-cell QC metrics to `adata.obs`

# %%
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)

# %% [markdown]
# ## 4. pre-filtering
#
# In this step, we apply a very conservative filter to the number of detected genes and the number of reads per cell. This step will remove cells and genes that are almost certainly non-informative, but still leaves many low-quality cells that need to be dealt with later.
#
# 1. remove cells that don't have a minimum number of counts and detected genes.

# %%
sc.pp.filter_cells(adata, min_counts=500)
sc.pp.filter_cells(adata, min_genes=200)

# %% [markdown]
# 2. remove genes that are not present in at least 20 cells.

# %%
sc.pp.filter_genes(adata, min_cells=20)

# %% [markdown]
# ## 5. outlier detection
#
# Similar to what is described in {cite}`germainPipeCompGeneralFramework2020` and {cite}`heumosBestPracticesSinglecell2023`,
# we apply a relativly lenient filter that removes cells that are outliers in at least two of the following categories
#  * log1p total counts by > 5 MAD
#  * log1p detected genes by > 5 MAD
#  * fraction of counts in top 20% of genes by > 5 MAD
#  * fraction of mitochondrial counts > 3 MAD
#
# :::{important}
# Using too stringent cutoffs may remove entire cell-types!
# It is better to apply permissive filtering here, and remove remaining low-quality cells during cell-type annotation.
# :::
#
# 1. Define outliers in each category

# %%
adata.obs["is_outlier_counts"] = aps.pp.is_outlier(adata, "log1p_total_counts", n_mads=5, groupby="sample")
adata.obs["is_outlier_genes"] = aps.pp.is_outlier(adata, "log1p_n_genes_by_counts", n_mads=5, groupby="sample")
adata.obs["is_outlier_top_20"] = aps.pp.is_outlier(adata, "pct_counts_in_top_20_genes", n_mads=5, groupby="sample")
adata.obs["is_outlier_mito"] = aps.pp.is_outlier(adata, "pct_counts_mt", n_mads=3, groupby="sample")

# %% [markdown]
# 2. Label cells that are outliers in at least two conditions in the `is_outlier` column.

# %%
adata.obs["is_outlier"] = (
    np.sum(
        adata.obs.loc[
            :,
            [
                "is_outlier_counts",
                "is_outlier_genes",
                "is_outlier_top_20",
                "is_outlier_mito",
            ],
        ],
        axis=1,
    )
    >= 2
)

# %% [markdown]
# ## 6. Subset data

# %%
adata.shape

# %%
adata_filtered = adata[~adata.obs["is_outlier"]].copy()

# %%
adata_filtered.shape

# %% [markdown]
# ## 7. Store result

# %%
adata_filtered.write_h5ad(f"{out_dir}/adata_filtered.h5ad")
