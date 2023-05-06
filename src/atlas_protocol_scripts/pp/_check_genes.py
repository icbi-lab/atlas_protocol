import re
from typing import List

import anndata
import numpy as np
import pandas as pd


def remove_gene_version(string):
    string = str(string)
    if string.startswith("ENSG"):
        return re.sub(r"\..*", "", string)
    else:
        return string


def append_duplicate_suffix(df: pd.DataFrame, column: str, sep: str) -> pd.DataFrame:
    """Appends a numeric suffix to each duplicated value in the specified column based on index position."""
    suffix_dict = {}
    df_sorted = df.sort_values(by=column)
    df_sorted = df_sorted.reset_index().sort_values(by="index").set_index("index")
    for i, val in enumerate(df_sorted[column]):
        if val in suffix_dict:
            suffix_dict[val] += 1
            df_sorted.at[df_sorted.index[i], column] = f"{val}{sep}{suffix_dict[val]}"
        else:
            suffix_dict[val] = 0
    return df_sorted.sort_index()


def find_unmapped_genes(adata: anndata.AnnData) -> List[str]:
    """Finds genes in the specified AnnData object that are not mapped to any ensembl id."""
    for column in ["Gene", "symbol"]:
        if column in adata.var.columns:
            unmapped = adata[:, adata.var_names == "nan"].var[column].tolist()
    return unmapped


def aggregate_duplicate_gene_ids(adata: anndata.AnnData, gene_names: List[str]) -> anndata.AnnData:
    """
    Collapse duplicate gene IDs in an AnnData object by summing their expression values.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to process.
    gene_names : List[str]
        A list of gene names to collapse by summing their expression values.

    Returns
    -------
    anndata.AnnData
        The modified AnnData object with collapsed gene IDs.

    Notes
    -----
    Gene IDs that are not in the input list will be kept in the output AnnData object.
    Gene IDs that are empty strings, "nan", or np.nan will be excluded from the input list.
    """
    # Remove empty values, "nan", and np.nan from the input gene_names list
    gene_names = [g for g in gene_names if g not in ["", "nan", np.nan]]
    if not gene_names:
        return adata

    # Subset for non-duplicate genes
    adata_tmp = adata[
        :, ~(adata.var_names.isin(gene_names))
    ].copy()  # not entirely sure that .isin() will work with symbols as substring will be matched, worst case the gene will appear as duplicate in the final adata
    # Subset for duplicate genes individually and compute sum of counts, create new AnnData objects with computed sums in .X
    sums = [np.sum(adata[:, gene].X, axis=1).A1 for gene in gene_names]
    adatas = [anndata.AnnData(X=sum.reshape((-1, 1)), obs=adata.obs.copy()) for sum in sums]

    # Concatenate the new AnnData objects and the original object along axis 1 to combine the computed sums with the original data
    for i, gene in enumerate(gene_names):
        adatas[i].var_names = [gene]
    adatas.append(adata_tmp)
    adata = anndata.concat(adatas, axis=1)
    return adata
