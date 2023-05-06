# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:CRCA-2023-crca-scanpy]
#     language: python
#     name: conda-env-CRCA-2023-crca-scanpy-py
# ---

# %% [markdown]
# # Merge datasets, harmonize annotations and metadata

# %% [markdown]
# Integrating single-cell RNA-seq datasets from multiple sources can provide numerous benefits, including increased statistical power, validation of findings across diverse conditions, and the identification of novel gene expression patterns that may be challenging to detect in individual datasets. However, the merging process presents two major challenges: harmonizing gene annotations and metadata across datasets to ensure consistency in downstream analyses.
#
# :::{note} gene annotations
#
# Ideally, access to raw FASTQ files would allow mapping to the same reference genome and annotations. However, in many cases, only processed data is available that may have been mapped to different genome annotations or versions. The two most commonly used gene annotation sources are GENCODE and Ensembl, which offer standardized gene models and annotations for various organisms. Best case scenario processed datasets have unique gene ids such as ensembl_ids, unfortunaetly often only gene symbols are provided that are not unique and can change across annotation versions and sources. sometimes provide only gene symbols.
# While it is possible to perform gene symbol-based integration, this approach is not always accurate, as gene symbols are not unique and can change between annotation versions. In contrast, before integrating the datasets we will map the available gene ids to the more consistent ensembl gene IDs that will enhance the accuracy and reproducibility of downstream analyses. :::
#
# Between different versions the ensembl gene ids will only change if the gene structure changes.
#
# Explain a bit more here? (e.g in newer versions new genes might be added, nothing we can do about it,
#                             if the gene id is the same the mapped gene region should have stayed the same. -> perfect!
#                             if the gene id has changed that means the gene structure has changed and we should not use it any more!)

# %% [markdown]
# ## 1. Load the required libaries and data

# %%
import anndata
import numpy as np
import pandas as pd
import scanpy as sc
import yaml
from scipy.sparse import csr_matrix

import atlas_protocol_scripts as aps

# %%
out_dir = "../../data/results/qc/"
# !mkdir -p {out_dir}

# %%
DATASETS = {
    "maynard_2020": "../../data/input_data_raw/maynard2020.h5ad",
    "lambrechts_2018": "../../data/input_data_raw/lambrechts_2018_luad_6653.h5ad",
    "ukim-v": "../../data/input_data_raw/ukim_v_batch1.h5ad",
}

# %%
datasets = {dataset_id: sc.read_h5ad(path) for dataset_id, path in DATASETS.items()}

# %%
# check that adata.X contains integers - requirement for scvi-tools integration
errors = {}
for name, adata in datasets.items():
    try:
        assert np.all(np.modf(adata.X.data)[0] == 0)
    except AssertionError:
        errors[name] = "X does not contain all integers"
errors

# %%
# round length corrected plate-based study
datasets["maynard_2020"].X.data = np.ceil(datasets["maynard_2020"].X.data).astype(int)

# %% [markdown]
# ## 2. Harmonize metadata

# %% [markdown]
# Before integrating the data we need to make sure to harmonize the metdata across our datasets. We will start by loading a pre-defined reference metadata yaml file that lists all columns we would like to have as well as the respective values that are allowed in every column. Using a helper function we can now quickly query what metadata is missing and if all values follow the same conventions.
#
# (Note: possible to use sfaira for metadata harmonization)

# %%
# Read the YAML file and load it into a dictionary
file_path = "../../tables/meta_reference.yaml"
with open(file_path) as f:
    ref_meta_dict = yaml.load(f, Loader=yaml.Loader)

# %%
# List reference columns from meta yaml file
ref_meta_cols = []
for key, _sub_dict in ref_meta_dict.items():
    ref_meta_cols.append(key)
ref_meta_cols

# %% [markdown]
# First we will check if all columns are present across all datasets

# %%
# loop over datasets and apply validate_obs function
invalid_columns = {}
for key, adata in datasets.items():
    try:
        aps.pp.validate_obs(adata, ref_meta_dict)
    except ValueError as e:
        invalid_columns[key] = e.args[0]
invalid_columns

# %%
# search reference dict for missing columns
aps.pp.search_dict(ref_meta_dict, ["condition"])

# %%
aps.pp.search_dict(ref_meta_dict, ["dataset", "platform"])

# %%
datasets["maynard_2020"].obs["dataset"] = "maynard_2020"
datasets["maynard_2020"].obs["platform"] = "smartseq2"

datasets["ukim-v"].obs["dataset"] = "ukim-v"
datasets["ukim-v"].obs["platform"] = "bd_rhapsody"
datasets["ukim-v"].obs["cell_type_salcher"] = "Unknown"

# %%
# loop over datasets and apply validate_obs function
invalid_columns = {}
for key, adata in datasets.items():
    try:
        aps.pp.validate_obs(
            adata,
            ref_meta_dict,
            keys_to_ignore=["dataset", "sample", "patient", "cell_type_salcher"],
        )
    except ValueError as e:
        invalid_columns[key] = e.args[0]
invalid_columns

# %% [markdown]
# We will ignore a few columns that contain unique ids as well as cell types - although we could of course define every allowed value in the yaml file.

# %%
# subset columns and keep only reference columns from meta yaml file
for adata in datasets:
    datasets[adata].obs = datasets[adata].obs[ref_meta_cols]

# %% [markdown]
# ## 3. Harmonize gene annotations

# %% [markdown]
# Before intgration we want ensembl ids without version numbers as var_names.
# note: we will have the best match between gene ids and symbols if we use the annotation that was used for mapping, can usually be found in the methods section of the paper or on GEO etc.

# %%
# load gtf for gene mapping
gtf_path = "../../tables/gencode.v32_gene_annotation_table.csv"
gtf = pd.read_csv(gtf_path)
gtf = aps.pp.append_duplicate_suffix(df=gtf, column="GeneSymbol", sep="-")
gene_ids = gtf.set_index("GeneSymbol")["Geneid"].to_dict()

# %%
datasets["lambrechts_2018"].var = datasets["lambrechts_2018"].var.rename_axis("symbol").reset_index()
datasets["lambrechts_2018"].var["ensembl"] = datasets["lambrechts_2018"].var["symbol"].map(gene_ids)
datasets["lambrechts_2018"].var_names = datasets["lambrechts_2018"].var["ensembl"].apply(aps.pp.remove_gene_version)

datasets["maynard_2020"].var.reset_index(inplace=True)
datasets["maynard_2020"].var_names = datasets["maynard_2020"].var["ensg"].apply(aps.pp.remove_gene_version)

datasets["ukim-v"].var.reset_index(inplace=True)
datasets["ukim-v"].var["ensembl"] = datasets["ukim-v"].var["Gene"].map(gene_ids)
datasets["ukim-v"].var_names = datasets["ukim-v"].var["ensembl"].apply(aps.pp.remove_gene_version)

# %%
# look how many genes were not mapped to ensembl ids
unmapped_dict = {}
for name, data in datasets.items():
    unmapped_genes = aps.pp.find_unmapped_genes(data)
    print(name, ":", len(unmapped_genes))
    unmapped_dict[name] = unmapped_genes

# %%
# remove genes without ensembl ids from the datasets
datasets["ukim-v"] = datasets["ukim-v"][:, ~(datasets["ukim-v"].var_names == "nan")]

# %%
# aggregate counts with the same id
for adata in datasets:
    duplicated_ids = datasets[adata].var_names[datasets[adata].var_names.duplicated()].unique()
    datasets[adata] = aps.pp.aggregate_duplicate_gene_ids(datasets[adata], duplicated_ids)
    assert datasets[adata].var_names.is_unique
    assert datasets[adata].obs_names.is_unique

# %%
# clean input data by removing not needed data
for col in ["counts_length_scaled", "tpm"]:
    del datasets["maynard_2020"].layers[col]

del datasets["ukim-v"].obsm["surface_protein"]

# %% [markdown]
# ## 4. Concat datasets to single adata

# %% [markdown]
# Finally the datasets are ready to be merged. We will also use the latest gene annotation from ensembl to update the gene ids and symbols. We could also use gencode.

# %%
# Outer join to keep all genes, fill_value=0 assuming that the removed gene expression was 0 or close to zero!
adata = anndata.concat(datasets, index_unique="_", join="outer", fill_value=0)

# %%
# Get latest ensembl annoation to update our genes
gtf_path = "../../tables/Homo_sapiens.GRCh38.109_gene_annotation_table.csv"
gtf = pd.read_csv(gtf_path)
gtf["ensembl"] = gtf["gene_id"].apply(aps.pp.remove_gene_version)
gtf["var_names"] = gtf["gene_name"].fillna(gtf["ensembl"])
gtf = aps.pp.append_duplicate_suffix(df=gtf, column="var_names", sep="-")

# %%
adata.var = pd.merge(
    pd.DataFrame({"ensembl": adata.var_names}),
    gtf,
    how="left",
    on="ensembl",
    validate="m:1",
).set_index("ensembl")

# Reorder by gtf (i.e. chromosome position)
gene_index = gtf[gtf["ensembl"].isin(adata.var_names)]["ensembl"].values
adata = adata[:, gene_index]

adata.var = adata.var.reset_index("ensembl")

adata.var_names = adata.var["var_names"].values
adata.var_names_make_unique()
del adata.var["var_names"]
adata.var_names.name = None

# %%
# Make sure samples are unique
adata.obs["sample"] = [f"{dataset}_{sample}" for dataset, sample in zip(adata.obs["dataset"], adata.obs["sample"])]

# %%
# Append dataset and sample info to barcodes
adata.obs_names = (
    adata.obs["dataset"].astype(str)
    + "_"
    + adata.obs["sample"].astype(str)
    + "_"
    + adata.obs_names.str.split("_").str[0]
)

# %%
assert adata.var_names.is_unique
assert adata.obs_names.is_unique

# %% [markdown]
# todo: have a threshold that removes genes if not present in at least 25 percent of studies?

# %% [markdown]
# ## 5. Store result

# %%
adata.write_h5ad(f"{out_dir}/adata.h5ad")
