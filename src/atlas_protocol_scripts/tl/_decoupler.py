import pandas as pd


def long_form_df_to_decoupler(
    df: pd.DataFrame,
    *,
    gene_col: str = "gene_id",
    p_col: str = "pvalue",
    log_fc_col: str = "log2FoldChange",
    group_col: str = "cell_type",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert a long-form data frame of DE results into a pair of logFC/pvalue matrices as used by decoupler.

    Parameters
    ----------
    df
        long-form pandas data frame
    gene_col
        column that contains the gene identifier (will be columns of the matrices)
    p_col
        column that contains the p-value or adjusted p-value
    log_fc_col
        column that contains the log fold change
    group_col
        column that contains some grouping information (will be rows of the matrices)
    """
    log_fc_mat = df.pivot(columns=gene_col, index=group_col, values=log_fc_col).fillna(0)
    p_mat = df.pivot(columns=gene_col, index=group_col, values=p_col).fillna(1)
    return log_fc_mat, p_mat
