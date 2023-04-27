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
# # Merge data and harmonized annotations

# %% [markdown]
# Integrating single-cell RNA-seq datasets from multiple sources can provide numerous benefits, including increased statistical power, validation of findings across diverse conditions, and the identification of novel gene expression patterns that may be challenging to detect in individual datasets. Ideally, access to raw FASTQ files would allow mapping to the same reference genome and annotations. However, in many cases, only processed data is available, which presents some challenges during data integration. For instance, the datasets may have been mapped to different genome annotations or versions, and sometimes provide only gene symbols.
# While it is possible to perform gene symbol-based integration, this approach is not always accurate, as gene symbols are not unique and can change between annotation versions. In contrast, before integrating the datasets we will map the available gene ids to the more consistent ensembl gene IDs that will enhance the accuracy and reproducibility of downstream analyses.
#
# The two most commonly used gene annotation sources are GENCODE and Ensembl, which offer standardized gene models and annotations for various organisms. Between different versions the ensembl gene ids will only change if the gene structure changes.
#
# Explain a bit more here? (e.g in newer versions new genes might be added, nothing we can do about it,
#                             if the gene id is the same the mapped gene region should have stayed the same. -> perfect!
#                             if the gene id has changed that means the gene structure has changed and we should not use it any more!)

# %% [markdown] tags=[]
# ## 1. Load the required libaries

# %%
# import atlas_protocol_scripts as aps
import anndata
import pandas as pd
import scanpy as sc

# %%
DATASETS = {
    "maynard_2020": "../../data/input_data_raw/maynard2020.h5ad",
    "lambrechts_2018": "../../data/input_data_raw/lambrechts_2018_luad_6653.h5ad",
    "ukim-v": "../../data/input_data_raw/ukim_v_batch1.h5ad",
}

# %%
datasets = {dataset_id: sc.read_h5ad(path) for dataset_id, path in DATASETS.items()}


# %%
# Get gene annotations from gtf file, remove Ensembl version number and append back sex chromosome info to new column “Ensembl”
def load_gtf(gtf_path):
    gtf = pd.read_csv(gtf_path, delimiter="\t", skipinitialspace=True, dtype={"Start": object, "End": object})
    gtf.insert(
        0,
        "Ensembl",
        gtf["Geneid"].str.replace(r"\.[^.]+$", "")
        + gtf["Geneid"].str.contains("_PAR_Y").replace({True: "_PAR_Y", False: ""}),
    )
    return gtf


# %%
gencode_v43 = load_gtf("../../tables/gencode.v43_gene_annotation_table.txt")
gencode_v33 = load_gtf("../../tables/gencode.v33_gene_annotation_table.txt")
gencode_v32 = load_gtf("../../tables/gencode.v32_gene_annotation_table.txt")


# %% [markdown]
# note: we will have the best match between gene ids and symbols if we use the annotation that was used for mapping, can usually be found in the methods section of the paper or on GEO etc.


# %%
def update_gene_annotation(adata, gtf, var_col="var_names", gtf_col="GeneSymbol"):
    """Annotates the genes in adata.var with their corresponding gene annotations from a gene annotation file
    and keeps the original columns in adata.var. Adata.var_names is set to “Ensembl” column from gtf file.

    Parameters
    ----------
    adata : AnnData object
    gtf : pandas.DataFrame object
        A data frame containing gene annotation information
    var_col : str, optional (default: “var_names”)
        The name of the column in adata.var to merge on
    gtf_col : str, optional (default: “GeneSymbol”)
        The name of the column in gtf to merge on

    Returns
    -------
    annotated AnnData object
        Contains the mapped genes and their corresponding annotations.

    Raises
    ------
    KeyError: If the “gtf_col” column is not present in the gene annotation data.
    ValueError: If the “var_col” column is not unique.
    """
    if gtf_col not in gtf.columns:
        raise KeyError(f"The gtf gene annotation data must contain a column named '{gtf_col}'")

    if var_col != "var_names":
        adata = adata.copy()
        var_index_name = adata.var_names.name
        adata.var["original_var_names"] = adata.var_names
        adata.var_names = adata.var[var_col]
    try:
        assert adata.var_names.is_unique
    except AssertionError:
        raise ValueError(f"'{var_col}' must be unique.")

    # Filter and remove duplicates from gene annotation data
    mapped_genes = gtf[gtf[gtf_col].isin(adata.var_names)].drop_duplicates(gtf_col, keep=False)
    mapped_genes = mapped_genes.reset_index(drop=True)

    # Filter adata for mapped/notmapped genes
    mapped_adata = adata[:, adata.var_names.isin(mapped_genes[gtf_col])].copy()

    mapped_adata.var = pd.merge(
        pd.DataFrame(mapped_adata.var).rename_axis("var_names").reset_index(),
        mapped_genes,
        how="left",
        left_on="var_names",
        right_on=[gtf_col],
        validate="m:1",
    ).set_index("Ensembl")

    if var_col != "var_names":
        del mapped_adata.var["var_names"]
        mapped_adata.var.rename(columns={"original_var_names": var_index_name}, inplace=True)

    return mapped_adata


# %%
params = [
    {"dataset": "maynard_2020", "gtf": gencode_v33, "var_col": "ensg", "gtf_col": "Geneid"},
    {"dataset": "lambrechts_2018", "gtf": gencode_v32, "var_col": "var_names", "gtf_col": "GeneSymbol"},
    {"dataset": "ukim-v", "gtf": gencode_v32, "var_col": "var_names", "gtf_col": "GeneSymbol"},
]

# %%
mapped_adata = {}
for param in params:
    dataset = param.pop("dataset")
    adata = datasets[dataset]
    mapped_adata[dataset] = update_gene_annotation(adata, **param)

# %% [markdown]
# # Have a look at not_mapped_genes + additional annotation of these

# %%
mapped = datasets["lambrechts_2018"].var_names.isin(mapped_adata["lambrechts_2018"].var["var_names"])
datasets["lambrechts_2018"][:, ~mapped].var

# %%
mapped = datasets["ukim-v"].var_names.isin(mapped_adata["ukim-v"].var["var_names"])
datasets["ukim-v"][:, ~mapped].var

# %% [markdown]
# todo: show how to add not mapped genes manually

# %% [markdown]
# # Concat datasets to single adata

# %%
# Outer join to keep all genes, fill_value=0 assuming that the removed gene expression was 0 or close to zero!
adata = anndata.concat(mapped_adata, join="outer", fill_value=0)

adata.var = pd.merge(
    pd.DataFrame({"Ensembl": adata.var_names}),
    gencode_v43,
    how="left",
    on="Ensembl",
    validate="m:1",
).set_index("Ensembl", drop=False)

adata.obs_names_make_unique()
assert adata.obs_names.is_unique

# Reorder by annotation and append back annotation info
assert adata.var_names.is_unique
adata = adata[:, gencode_v43.loc[gencode_v43["Ensembl"].isin(adata.var_names), "Ensembl"].values]

adata.var_names = adata.var["GeneSymbol"]
adata.var_names_make_unique()
adata.var_names.name = None

# %% [markdown]
# -> now we could for example only keep genes that are protein coding

# %%
adata_protein = adata[:, adata.var["Class"] == "protein_coding"].copy()

# %%
adata_protein.var

# %% [markdown]
# todo: remove adata.obs columns not needed in downstream analysis

# %% [markdown]
# todo: have a threshold that removes genes if not present in at least 25 percent of studies?
