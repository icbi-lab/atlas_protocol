import re

from anndata import AnnData
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix


def remove_gene_version(string: str) -> str:
    """Remove gene version from ensembl gene id that start with "ENSG"."""
    string = str(string)
    if string.startswith("ENSG"):
        return re.sub(r"\..*", "", string)
    else:
        return string


def append_duplicate_suffix(df: pd.DataFrame, column: str, sep: str) -> pd.DataFrame:
    """Appends a numeric suffix to each duplicated value in the specified column based on index position."""
    suffix_dict = {}
    df_sorted = df.sort_values(by=column).sort_index()
    for i, val in enumerate(df_sorted[column]):
        if val in suffix_dict:
            suffix_dict[val] += 1
            df_sorted.at[df_sorted.index[i], column] = f"{val}{sep}{suffix_dict[val]}"
        else:
            suffix_dict[val] = 0
    return df_sorted.sort_index()


def find_unmapped_genes(adata: AnnData, column: str = "var_names") -> list[str]:
    """Finds genes in the specified AnnData object that are not mapped to any ensembl id."""
    unmapped = []
    if column in adata.var.columns:
        col_data = adata.var[column]
    else:
        col_data = adata.var_names
    unmapped += col_data[~adata.var_names.str.startswith("ENSG")].tolist()
    return unmapped


def aggregate_duplicate_gene_ids(adata: AnnData, gene_names: list[str]) -> AnnData:
    """Collapse duplicate gene IDs in an AnnData object by summing their expression values.

    Parameters
    ----------
    adata
        AnnData object
    gene_names
        list of gene names to collapse by summing their expression values.
    """
    gene_names = [g for g in gene_names if g not in ["", "nan", np.nan]]
    if not gene_names:
        return adata

    adata.X = adata.X.tocsc()
    sums = [csc_matrix.sum(adata[:, gene].X, axis=1) for gene in gene_names]

    # Reshape the list of arrays to match the dimensions of adata.X
    new_X = np.concatenate(sums).reshape(-1, len(gene_names))

    adata_sums = anndata.AnnData(
        X=new_X,
        obs=adata.obs,
        var=pd.DataFrame(gene_names).set_index(0).rename_axis(None),
    )
    # Remove duplicated genes before concat
    adata = adata[:, ~(adata.var_names.isin(gene_names))]
    adata_tmp = anndata.concat([adata, adata_sums], axis=1)
    adata_tmp.obs = adata.obs  # Hope this preserves the right order ...
    return adata_tmp
