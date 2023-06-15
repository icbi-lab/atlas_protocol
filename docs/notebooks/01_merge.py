# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python [conda env:CRCA-2023-crca-scanpy]
#     language: python
#     name: conda-env-CRCA-2023-crca-scanpy-py
# ---

# %% [markdown]
# # Merge datasets, harmonize annotations and metadata
#
# Integrating single-cell RNA-seq datasets from multiple sources can provide numerous benefits, including increased statistical power, validation of findings across diverse conditions, and the identification of novel gene expression patterns that may be challenging to detect in individual datasets. However, the merging process presents two major challenges: harmonizing gene annotations and metadata across datasets to ensure consistency in downstream analyses.
#
#
# :::{Important}
# **types of metadata**
#
# There are essentially four different levels of metadata that need to be gathered and organized:
#
#  * study metadata
#  * patient metadata
#  * sample metadata
#  * cell metadata
#
# Commonly, this information can be found in the original publication, either in the methods section or supplementary tables. Some metadata might also be available on data repositories such as [GEO](https://www.ncbi.nlm.nih.gov/geo/), [SRA](https://www.ncbi.nlm.nih.gov/sra), [Zenodo](https://zenodo.org), [cellxgene](https://cellxgene.cziscience.com), [synapse](https://www.synapse.org), or custom webpages that were created for a specific study. As a last resort, you can contact the authors for assistance.
# :::
#
# :::{See also}
# - An example of how to organize metadata across multiple datasets can be found in the study conducted by {cite}`zilbauer2023`.
# - [sfaira](https://sfaira.readthedocs.io/en/latest/): Framework to organize and use publicly available datasets: {cite}`fischer2021sfaira`
# :::
#
# :::{Study metadata}
# Typical examples:
#
#  - sequencing platform: Smartseq vs 10x 3' v1/v2/v3, 10x 5', etc.
#  - tissue prcessing: fresh vs. frozen
#  - protocol: single cell vs. single nuclei
#  - prior cell type enrichment: naive, CD45+, etc.
#
# In addition, it is worth noting whether there are supplementary assays, such as VDJ-seq and/or CITE-seq data, that can help resolve immune cell heterogeneity
# :::
#
# :::{Patient metadata}
# Typical examples:
#
#  - sex: male vs. female
#  - age
#  - ethnicity
#  - treatment status: naive vs. treated
# :::
#
# :::{Sample metadata}
# Typical examples:
#
#  - sample type: tumor/normal/pbmc/metastasis/etc.
#  - sample tissue: colon/liver/lymph node/blood/etc.
#  - primary tumor location
#  - tumor stage: e.g. TNM system
#  - histological tumor type
#  - known driver mutations
# :::
#
# :::{Cell metadata}
# If available, the cell type annotation from the original study can be used for reference mapping. See <!-- [link to chapter](TODO) -->
# :::

# %% [markdown]
# ## 1. Load the required libaries and datasets

# %%
import anndata
import atlas_protocol_scripts as aps
import numpy as np
import pandas as pd
import scanpy as sc
import yaml

# %%
out_dir = "../../data/results/merge/"
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
# Check that adata.X contains integers - requirement for scvi-tools integration
errors = {}
for name, adata in datasets.items():
    try:
        assert np.all(np.modf(adata.X.data)[0] == 0)
    except AssertionError:
        errors[name] = "X does not contain all integers"
errors

# %%
# Round length corrected plate-based study
datasets["maynard_2020"].X.data = np.ceil(datasets["maynard_2020"].X.data).astype(int)

# %% [markdown]
# ## 2. Harmonize metadata
#
# To ensure streamlined metadata across our datasets, we will use a custom reference metadata YAML file that specifies the desired columns and the permissible values for each column in `adata.obs`. Here is a shortened example of different key-value pairs along with descriptions:
#
# ```{yaml}
# origin:
#     values:
#         - tumor_primary
#         - normal_adjacent
#         - tumor_edge
#         - tumor_middle
#         - tumor_metastasis
#         - nan
#     description: Sample origin
# condition:
#     values:
#         - LUAD
#         - LSCC
#         - NSCLC
#     description:
#         Lung adenocarcinoma (LUAD) and lung squamous cell carcinoma (LSCC)
#         are the most common subtypes of non-small-cell lung cancer (NSCLC)
# platform:
#     values:
#         - 10x_3p_v2
#         - smartseq2
#         - bd_rhapsody
#     description: Protocol that was used for single cell sequencing
# ```
#
# The reference metadata YAML file serves as the primary location to define key-value pairs for different metadata columns, and any additional metadata columns can be easily added. It enables easy querying of allowed values during metadata collection. Furthermore, we will use it as a final check to ensure that all columns in the merged `adata.obs` follow the defined conventions using a helper function {func}`~atlas_protocol_scripts.pp.validate_obs`.

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

# %%
# Loop over datasets and apply validate_obs function to check if all columns are present across all datasets
for key, adata in datasets.items():
    try:
        aps.pp.validate_obs(adata.obs, ref_meta_dict)
    except ValueError as e:
        raise ValueError(e.args[0])

# %% [markdown]
# The `ValueError` tells us that we need to add missing metadata columns in some of the datasets.

# %% tags=[]
# Search reference dict for permissible values of missing columns
ref_meta_dict["platform"]

# %%
# Add missing metadata: we will need "dataset" to make patient and sample ids unique; "platform" to check how well the integration worked; "cell_type_salcher" for seed annotation.
datasets["maynard_2020"].obs["dataset"] = "maynard_2020"
datasets["maynard_2020"].obs["platform"] = "smartseq2"

datasets["ukim-v"].obs["dataset"] = "ukim-v"
datasets["ukim-v"].obs["platform"] = "bd_rhapsody"
datasets["ukim-v"].obs["cell_type_salcher"] = "Unknown"

# %%
# Loop over datasets and apply validate_obs function. Additionally, we will exclude columns from the permissible values check that are expected to be unique within each dataset.
for key, adata in datasets.items():
    try:
        aps.pp.validate_obs(
            adata.obs,
            ref_meta_dict,
            keys_to_ignore=["dataset", "sample", "patient", "cell_type_salcher"],
        )
    except ValueError as e:
        raise ValueError(e.args[0])

# %%
# Subset adata.obs columns to keep only reference columns from meta yaml file
for adata in datasets:
    datasets[adata].obs = datasets[adata].obs[ref_meta_cols]

# %% [markdown]
# ## 3. Harmonize gene annotations
#
# Ideally, access to raw FASTQ files would allow mapping to the same reference genome and annotations. However, in many cases, only processed data is available that may have been mapped to different genome annotations or versions. The two most commonly used gene annotation sources are [GENCODE](https://www.gencodegenes.org) and [Ensembl](https://www.ensembl.org/index.html), which offer standardized gene models and annotations for various organisms.
#
# While it is possible to perform gene symbol-based integration, this approach is not always accurate, as gene symbols are not unique and can change between annotation versions. In contrast, before integrating the datasets we will map the available gene ids to the more consistent ensembl gene IDs that will enhance the accuracy and reproducibility of downstream analyses.
#
# :::{note}
# **ENSG id conversion between versions:**
#
#  - If the ENSG ID remains unchanged, it signifies that the gene structure remains the same, while updated structures will have different ENSG IDs.
#  - Newer annotation versions often include additional annotated transcripts. It's important to note that most of these transcripts are non-protein coding. Without access to fastq files, it is not possible to update the annotations to incorporate these newly annotated transcripts.
#  - Deprecated ENSG IDs will be absent from newer annotations.
#
# Once merged, we will update all annotations using the most recent reference, and any deprecated IDs will be removed to obtain the best possible updated version.
# :::
#
# :::{See also}
#  - {cite}`bruford2020`
# :::
#
# :::{Alternative approaches}
#  - [python gtfparse](https://github.com/openvax/gtfparse) is an alternative tool for loading a gtf file into Python
#  - When using gene symbols for mapping, it is recommended to update them using libraries such as [mygene-py](https://docs.mygene.info/projects/mygene-py/en/latest/) or by using the ensembl biomart aliases [ensembl biomart aliases](https://biomart.genenames.org) before merging the datasets.
# :::

# %%
# Load reference gtf for gene mapping, and create a dictonary with symbol-ENSG ID pairs.
gtf_path = "../../tables/gencode.v32_gene_annotation_table.csv"
gtf = pd.read_csv(gtf_path)

# When making var_names unique, both Scanpy and Seurat typically append a sequential number to duplicate gene symbols.
# To match all symbol-ENSG ID pairs, we need to emulate this sequential numbering approach.
gtf = aps.pp.append_duplicate_suffix(df=gtf, column="GeneSymbol", sep="-")
gene_ids = gtf.set_index("GeneSymbol")["Geneid"].to_dict()

# %%
# 1. Map the dictonary to available symbol annotation and fill missing keys with the respective symbol to create new column "ensembl".
# 2. Remove Ensembl ID version numbers and set column "ensembl" as var_names.

datasets["lambrechts_2018"].var = datasets["lambrechts_2018"].var.rename_axis("symbol").reset_index()
datasets["lambrechts_2018"].var["ensembl"] = (
    datasets["lambrechts_2018"].var["symbol"].map(gene_ids).fillna(value=datasets["lambrechts_2018"].var["symbol"])
)
datasets["lambrechts_2018"].var_names = datasets["lambrechts_2018"].var["ensembl"].apply(aps.pp.remove_gene_version)

datasets["maynard_2020"].var.reset_index(inplace=True)
datasets["maynard_2020"].var_names = datasets["maynard_2020"].var["ensg"].apply(aps.pp.remove_gene_version)

datasets["ukim-v"].var.reset_index(inplace=True)
datasets["ukim-v"].var["ensembl"] = (
    datasets["ukim-v"].var["Gene"].map(gene_ids).fillna(value=datasets["ukim-v"].var["Gene"])
)
datasets["ukim-v"].var_names = datasets["ukim-v"].var["ensembl"].apply(aps.pp.remove_gene_version)

# %%
# Look how many genes were not mapped to ensembl ids
unmapped_dict = {}
for name, data in datasets.items():
    unmapped_genes = aps.pp.find_unmapped_genes(data)
    print(name, ":", len(unmapped_genes))
    unmapped_dict[name] = unmapped_genes

# %%
# Remove genes without ensembl ids from the datasets
datasets["ukim-v"] = datasets["ukim-v"][:, (datasets["ukim-v"].var_names.str.startswith("ENSG"))]

# %% [markdown]
# :::{note}
# To achieve the best match between ENSG IDs and gene symbols, it is advisable to use the annotation that was originally used for mapping. This information is typically available in the methods section of the paper or can be obtained from the associated data repository. If it is unavailable, an alternative approach is to deduce the annotation by downloading different versions and checking the number of unmapped genes after mapping.
# :::

# %%
# Aggregate counts with the same id
for adata in datasets:
    duplicated_ids = datasets[adata].var_names[datasets[adata].var_names.duplicated()].unique()
    datasets[adata] = aps.pp.aggregate_duplicate_gene_ids(datasets[adata], duplicated_ids)
    assert datasets[adata].var_names.is_unique
    assert datasets[adata].obs_names.is_unique

# %%
# Clean input data by removing not needed data
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
# Make sure samples are unique
adata.obs["sample"] = [f"{dataset}_{sample}" for dataset, sample in zip(adata.obs["dataset"], adata.obs["sample"])]

# Append dataset and sample info to barcodes
adata.obs_names = (
    adata.obs["dataset"].astype(str)
    + "_"
    + adata.obs["sample"].astype(str)
    + "_"
    + adata.obs_names.str.split("_").str[0]
)

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

# Put up to date gene symbols as var_names
adata.var_names = adata.var["var_names"].values
adata.var_names_make_unique()
del adata.var["var_names"]
adata.var_names.name = None

# %%
# Filter genes: must be expressed in at least 25 percent of the samples
adata.obs["sample"] = pd.Categorical(adata.obs["sample"], categories=adata.obs["sample"].unique())

res = pd.DataFrame(columns=adata.var_names, index=adata.obs["sample"].cat.categories)
for sample in adata.obs["sample"].cat.categories:
    res.loc[sample] = adata[adata.obs["sample"].isin([sample]), :].X.sum(0)

keep = res.columns[res[res == 0].count(axis=0) / len(res.index) >= 0.25]

# %% [markdown]
# :::{important}
# Check your gene filter cut-off, particularly when dealing with a limited number of samples/datasets and/or studies with prior cell type enrichment. Setting an overly stringent filter may result in the loss of important marker genes that could be valuable for downstream analyses.
# :::

# %%
# Subset adata to remove genes that dont pass the cut-off
adata = adata[:, keep].copy()

# %%
assert adata.var_names.is_unique
assert adata.obs_names.is_unique

# %%
# Look at final adata
adata

# %% [markdown]
# ## 5. Store result

# %%
adata.write_h5ad(f"{out_dir}/adata.h5ad")
