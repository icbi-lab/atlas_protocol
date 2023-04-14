import warnings
from functools import reduce
from operator import and_
from typing import Sequence, Union

import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData, ImplicitModificationWarning


def pseudobulk(adata, *, groupby: Union[str, Sequence[str]], aggr_fun=np.sum, min_obs=10, log_norm=False) -> AnnData:
    """Calculate Pseudobulk of groups.

    Parameters
    ----------
    adata
        annotated data matrix
    groupby
        One or multiple columns to group by
    aggr_fun
        Callback function to calculate pseudobulk. Must be a numpy ufunc supporting
        the `axis` attribute.
    min_obs
        Exclude groups with less than `min_obs` observations
    log_norm
        If True, compute log(CPM) of the result. Use this only when raw counts are in adata.X

    Returns
    -------
    New anndata object with same vars as input, but reduced number of obs.
    """
    if isinstance(groupby, str):
        groupby = [groupby]

    combinations = adata.obs.loc[:, groupby].drop_duplicates()

    if adata.is_view:
        # for whatever reason, the pseudobulk function is terribly slow when operating on a view.
        adata = adata.copy()

    # precompute masks
    masks = {}
    for col in groupby:
        masks[col] = {}
        for val in combinations[col].unique():
            masks[col][val] = adata.obs[col] == val

    expr_agg = []
    obs = []

    for comb in combinations.itertuples(index=False):
        mask = reduce(and_, (masks[col][val] for col, val in zip(groupby, comb)))
        if np.sum(mask) < min_obs:
            continue
        expr_row = aggr_fun(adata.X[mask, :], axis=0)
        obs_row = comb._asdict()
        obs_row["n_obs"] = np.sum(mask)
        # convert matrix to array if required (happens when aggregating spares matrix)
        try:
            expr_row = expr_row.A1
        except AttributeError:
            pass
        obs.append(obs_row)
        expr_agg.append(expr_row)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ImplicitModificationWarning)
        adata = AnnData(
            X=np.vstack(expr_agg),
            var=adata.var,
            obs=pd.DataFrame.from_records(obs),
        )

        if log_norm:
            sc.pp.normalize_total(adata, target_sum=1e6)
            sc.pp.log1p(adata)

        return adata
