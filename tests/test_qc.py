import pandas as pd
from anndata import AnnData

import atlas_protocol_scripts as aps


def test_package_has_version():
    assert aps.__version__ is not None


def test_is_outlier():
    adata = AnnData(
        obs=pd.DataFrame().assign(
            x=[0, 1, 20, 20, 21, 21, 22, 22, 40, 41], group=["a", "a", "b", "b", "b", "b", "b", "b", "a", "a"]
        )
    )
    assert aps.pp.is_outlier(adata, "x", groupby=None).tolist() == [
        True,
        True,
        False,
        False,
        False,
        False,
        False,
        False,
        True,
        True,
    ]
    assert aps.pp.is_outlier(adata, "x", groupby=None, n_mads=100).tolist() == [False] * 10
    assert aps.pp.is_outlier(adata, "x", groupby="group", n_mads=0.5).tolist() == [
        True,
        True,
        True,
        True,
        False,
        False,
        True,
        True,
        True,
        True,
    ]
