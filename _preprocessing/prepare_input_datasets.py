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
# Prepare input datasets for the tutorial.
# We choose datasets from three different platforms:
#  * Lambrechts_Thienpont_2018_6653 (10x 3' v2)
#  * Maynard_Bivona_2020 (SS2)
#  * UKIM-V (BD)
#  * UKIM-V-2 (BD, for scarches)

# %%
import scanpy as sc

# %%
atlas = sc.read_h5ad(
    "/home/sturm/projects/2020/pircher-scrnaseq-lung/data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad"
)

# %% [markdown]
# ### Copy UKIM-V dataset without change

# %%
ukimv1 = sc.read_h5ad("../data/input_data_raw/ukim_v_batch1.h5ad")

# %%
ukimv2 = sc.read_h5ad("../data/input_data_raw/ukim_v_batch2.h5ad")

# %%
ukimv1.write_h5ad("../data/input_data_zenodo/ukim_v_batch1.h5ad")

# %%
ukimv2.write_h5ad("../data/input_data_zenodo/ukim_v_batch2.h5ad")

# %% [markdown]
# ### Add known cell-type annotations to Maynard/Lambrechts (to skip seed annotation)

# %% [markdown]
# #### Maynard

# %%
maynard = sc.read_h5ad("../data/input_data_raw/maynard2020.h5ad")

# %%
obs_maynard = atlas[atlas.obs["dataset"] == "Maynard_Bivona_2020"].obs
obs_maynard.index = obs_maynard.index.str.replace("-\\d+$", "", regex=True)

# %%
maynard.obs["cell_type_salcher"] = obs_maynard["cell_type"]

# %%
maynard.write_h5ad("../data/input_data_zenodo/maynard2020.h5ad")

# %%
# sc.pp.highly_variable_genes(obs_maynard, n_top_genes=2000, flavor="seurat_v3")
# sc.pp.normalize_total(obs_maynard)
# sc.pp.log1p(obs_maynard)

# %% [markdown]
# #### Lambrechts

# %%
lambrechts = sc.read_h5ad("../data/input_data_raw/lambrechts_2018_luad_6653.h5ad")

# %%
obs_lambrechts = atlas[atlas.obs["dataset"] == "Lambrechts_Thienpont_2018_6653"].obs
obs_lambrechts.index = obs_lambrechts.index.str.replace("-\\d+$", "", regex=True)

# %%
lambrechts.obs["cell_type_salcher"] = obs_lambrechts["cell_type"]

# %%
lambrechts.write_h5ad("../data/input_data_zenodo/lambrechts_2018_luad_6653.h5ad")

# %%
# sc.pp.highly_variable_genes(lambrechts, n_top_genes=2000, flavor="seurat_v3")
# sc.pp.normalize_total(lambrechts)
# sc.pp.log1p(lambrechts)
# sc.tl.pca(lambrechts)
# sc.pp.neighbors(lambrechts)
# sc.tl.umap(lambrechts)

# %%
# sc.pl.umap(lambrechts, color=["patient", "cell_type"])

# %%
# sc.tl.leiden(lambrechts)

# %%
# sc.pl.umap(lambrechts, color="leiden")
