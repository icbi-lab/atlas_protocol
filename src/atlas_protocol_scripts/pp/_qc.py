from typing import Optional

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy.stats import median_abs_deviation


def is_outlier(adata: AnnData, metric_col: str, *, groupby: Optional[str] = None, n_mads: float = 5) -> pd.Series:
    """Detect outliers by median absolute deviation (MAD).

    Adapted from https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#motivation

    Parameters
    ----------
    adata
        AnnData object
    metric_col
        column in adata.obs to consider
    groupby
        grouping variable. If specified, outliers will be determined for each value of this grouping variable separately.
        E.g. `dataset`.
    n_mads
        label points that are outside of median +/- nmads * MAD.
    """

    def _is_outlier(df):
        """Check if metric value deviates from the median by more than n_mads MADs."""
        metric_values = df[metric_col]
        return np.abs(metric_values - np.median(metric_values)) > n_mads * median_abs_deviation(metric_values)

    if groupby is None:
        return _is_outlier(adata.obs)
    else:
        return adata.obs.groupby(groupby).apply(_is_outlier).droplevel(groupby).reindex(adata.obs_names)  # type: ignore
