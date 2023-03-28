---
jupyter:
    jupytext:
        text_representation:
            extension: .md
            format_name: markdown
            format_version: "1.3"
            jupytext_version: 1.14.4
    kernelspec:
        display_name: Python [conda env:.conda-pircher-sc-integrate2]
        language: python
        name: conda-env-.conda-pircher-sc-integrate2-py
---

Create demo dataset from full atlas including the three datasets in our example that can be used for downstream analyses

```python
import scanpy as sc
```

```python
adata = sc.read_h5ad(
    "/home/sturm/projects/2020/pircher-scrnaseq-lung/data/20_build_atlas/add_additional_datasets/03_update_annotation/artifacts/full_atlas_merged.h5ad"
)
```

```python
adata.obs["dataset"].value_counts()
```

```python
adata_sub = adata[
    adata.obs["dataset"].isin(
        ["Maynard_Bivona_2020", "Lambrechts_Thienpont_2018_6653", "UKIM-V"]
    ),
    :,
].copy()
```

```python
adata_sub
```

```python
sc.pp.neighbors(adata_sub, use_rep="X_scANVI")
```

```python
sc.tl.umap(adata_sub, init_pos="X_umap")
```

```python
sc.pl.umap(adata_sub, color="cell_type")
```

```python
adata_sub.write_h5ad("../data/input_data_zenodo/atlas-integrated-annotated.h5ad")
```

```python

```
