# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:.conda-pircher-sc-integrate2]
#     language: python
#     name: conda-env-.conda-pircher-sc-integrate2-py
# ---

# %% [markdown]
# # Compositional analysis

# %%
import altair as alt
import scanpy as sc

# %%
input_adata = "../../data/input_data_zenodo/atlas-integrated-annotated.h5ad"

# %%
adata = sc.read_h5ad(input_adata)

# %%
sc.pl.umap(adata, color="cell_type_coarse")

# %%
cell_type_fractions = (
    adata.obs.loc[(adata.obs["origin"] == "tumor_primary") & adata.obs["condition"].isin(["LUAD", "LUSC"]), :]
    .groupby(["condition", "patient", "cell_type_coarse"])
    .size()
    .reset_index(name="n")
    .groupby(["condition", "patient"])
    .apply(lambda x: x.assign(frac=x["n"] / x["n"].sum()))
    .dropna()
    .groupby(["condition", "cell_type_coarse"])
    .agg("mean")
    .reset_index()
)

# %%
alt.Chart(cell_type_fractions).encode(x="condition", y="frac", color="cell_type_coarse").mark_bar()

# %%
