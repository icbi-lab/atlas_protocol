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
# # Cell-type annotation of "seed" datasets
#
# scVI and scANVI are variational autoencoders that embed the single-cell expression data into a low-dimensional latent space, while removing batch effects. This is what we will be doing in chapter {ref}`data_integration`. While scVI is an unsupervised method that only considers the count data, scANVI is a "semi-supervised" data that takes into account known cell-type labels of one or multiple datasets.
#
# In an independent benchmark, the semi-supervised variant scANVI has outperformed scVI and other methods for atlas-level data integration {cite}`lueckenBenchmarkingAtlaslevelData2022`.
#
# In order to leverage the scANVI algorithm for building the atlas, we are going to prepare cell-type labels for two datasets with very different characteristics:
#  * `Lambrechts_Thienpont_2018_6653`, which has been sequenced using the dropblet-based, UMI-corrected 10x Genomics 3' v2 protocol
#  * `Maynard_Bivona_2020`, which has been sequenced using the well-based full-length Smart-seq2 protocol {cite}`picelliSmartseq2SensitiveFulllength2013`.

# %% [markdown]
# ## 1. Import the required libraries

# %%
import scanpy as sc

# %% [markdown]
# ## 2. Load the input data

# %%
lambrechts2018 = sc.read_h5ad("../../data/input_data_zenodo/lambrechts_2018_luad_6653.h5ad")
maynard2020 = sc.read_h5ad("../../data/input_data_zenodo/maynard2020.h5ad")

# %% [markdown]
# ## 3. Define and create output directory

# %%
out_dir = "../../results/seed_annotation"

# %%
# !mkdir -p {out_dir}

# %% [markdown]
# ## 3. preprocess each dataset individually
#
# TODO either based on scVI or just using normalize_total/log1p. Do this once filtering is complete.

# %% [markdown]
# ## 4. annotate cell-types for each dataset individually
#
# Seed datasets can be annotated based on unsupervised clustering and marker genes as shown in section {ref}`cell_type_annotation`. For the sake of this tutorial, we simply re-use the cell-type annotations from {cite}`salcherHighresolutionSinglecellAtlas2022a`.

# %%
lambrechts2018.obs["cell_type"] = lambrechts2018.obs["cell_type_salcher"]
maynard2020.obs["cell_type"] = maynard2020.obs["cell_type_salcher"]

# %% [markdown]
# ## 5. Store annotated AnnData objects

# %%
lambrechts2018.write_h5ad(f"{out_dir}/lambrechts_2018_luad_6653_annotated.h5ad")
maynard2020.write_h5ad(f"{out_dir}/maynard2020_annotated.h5ad")
