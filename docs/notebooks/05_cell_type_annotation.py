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
# (cell_type_annotation)=
# # Cell-type annotation
#
# Cell-type annotation is fundamental for data interpretation and a prerequisite for almost any
# downstream analysis. Performing cell-type annotation manually based on clustering is a widely-used
# and straightforward approach. To this end, previously known marker genes can be checked
# for their expression across clusters and cell-type labels assigned accordingly. Alternatively, as a
# more unbiased approach, it is possible to perform differential gene expression analysis between clusters, and assign cell-types based on cluster-specific differentially expressed
# (DE) genes. Panels of marker genes can be obtained from publications or cell-type databases such as CellMarker
# [103] or PanglaoDB [104].

# %% [markdown]
# ## 1. Load the required libraries

# %%
import pandas as pd
import scanpy as sc

# %% [markdown]
# ## 2. Load input data

# %%
adata = sc.read_h5ad("../../data/input_data_zenodo/atlas-integrated-annotated.h5ad")

# %%
markers = pd.read_csv("../../data/input_data_zenodo/cell_type_markers_lung.tsv", sep="\t")

# %%
markers

# %% [markdown]
# ## 3. Perform unsupervised clustering using the leiden algorithm

# %%
sc.tl.leiden(adata, resolution=1)

# %%
sc.pl.umap(adata, color="leiden")

# %%
