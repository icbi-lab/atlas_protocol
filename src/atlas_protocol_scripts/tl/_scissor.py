import warnings

import numpy as np
import pandas as pd
import scipy.stats


def scissor_wilcoxon_test(df: pd.DataFrame, pos_col: str = "scissor+", neg_col: str = "scissor-") -> pd.Series:
    """Tests if the fractions of scissor+ cells is significantly different from the fraction of scissor- cells.

    Applies a wilcoxon test on the scissor+ vs. scissor- fractions.

    Parameters
    ----------
    df
        a data frame where each row is a biological replicate (e.g. patient) and there are at least
        two columns with the fractions of scissor+ and scissor- cells, respectively.
    pos_col
        column in df that contains the fraction of scissor+ cells
    neg_col
        column in df that contains the fraction of scissor- cells

    Returns
    -------
    a Series with the following columns:
        * `pos_col`: the mean fraction across all patients in `pos_col`
        * `neg_col`: the mean fraction across all patisn in `neg_col`
        * pvalue: the p-value as computed using the wilcoxon test
        * `log2_ratio`: The log2 ratio of mean fractions. Computed as `log2(mean(pos_col)) - log2(mean(neg_col))`.
    """
    with warnings.catch_warnings():
        # Let's ignore: UserWarning: Exact p-value calculation does not work if there are zeros. Switching to normal approximation
        warnings.simplefilter("ignore")
        _, p = scipy.stats.wilcoxon(
            df[pos_col].values,
            df[neg_col].values,
            # deal with zeros
            zero_method="zsplit",
        )
    with np.errstate(divide="ignore"):
        return pd.Series(
            {
                pos_col: df[pos_col].mean(),
                neg_col: df[neg_col].mean(),
                "pvalue": p,
                "log2_ratio": np.log2(df[pos_col].mean()) - np.log2(df[neg_col].mean()),
            }
        )
