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

# %% tags=["hide-cell"]
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# # Quality control and filtering
#
# ```{note}
# **datasets vs. studies**: By “study”, we refer to a scientific publication, while with “dataset”,
# we refer to a set of samples that was generated using the same sequencing platform and processed in the same way.
# One study may contain one or multiple datasets. Datasets must be processed independently.
# ```

# %%
import matplotlib.pyplot as plt
import scanpy as sc
from IPython.core.display import HTML, display

import atlas_protocol_scripts as aps

# %%
DATASETS = {
    "maynard_2020": "../../data/input_data_raw/maynard2020.h5ad",
    "lambrechts_2018": "../../data/input_data_raw/lambrechts_2018_luad_6653.h5ad",
    "ukim-v": "../../data/input_data_raw/ukim_v_batch1.h5ad",
}
THRESHOLDS = {
    "maynard_2020": {
        "min_genes": 600,
        "min_counts": 20000,
        "max_counts": 20000000,
        "max_pct_mito": 30,
    },
    "lambrechts_2018": {
        "min_genes": 250,
        "min_counts": 1200,
        "max_counts": 40000,
        "max_pct_mito": 20,
    },
    "ukim-v": {
        "min_genes": 200,
        "min_counts": 2000,
        "max_counts": 100000,
        "max_pct_mito": 30,
    },
}

# %%
datasets = {dataset_id: sc.read_h5ad(path) for dataset_id, path in DATASETS.items()}

# %%
for dataset, thresholds in THRESHOLDS.items():
    display(HTML(f"<h2>{dataset}</h2>"))
    adata = datasets[dataset]
    if "mito" not in adata.var.columns:
        adata.var["mito"] = adata.var_names.str.lower().str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=("mito",), log1p=False, inplace=True, percent_top=None)
    figwidth = min(max(adata.obs["sample"].unique().size * 0.5, 2), 20)
    fig, ax = plt.subplots(1, 1, figsize=(figwidth, 5))
    sc.pl.violin(adata, "total_counts", groupby="sample", rotation=90, log=True, cut=0, ax=ax)
    aps.pl.plot_qc_metrics(adata, **thresholds)
    aps.pl.plot_qc_metrics(adata, cumulative=True, **thresholds)

# %%
