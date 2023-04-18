# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:.conda-atlas_protocol_23]
#     language: python
#     name: conda-env-.conda-atlas_protocol_23-py
# ---

# %% [markdown]
# # Differential gene expression analysis per cell-type between conditions

# %%
# Load libraries
import glob
import os
import subprocess as sp
import warnings

# Libraries for visualization
from typing import Sequence

import decoupler as dc
import numpy as np
import pandas as pd
import scanpy as sc

warnings.filterwarnings("ignore")
from itertools import zip_longest
from math import ceil

import matplotlib.pyplot as plt
import scipy.stats
import seaborn as sns
import statsmodels.stats.multitest
from IPython.display import display

# set PATH env variable to conda env for specific R version.
# To use [DESeq2, R version "4.2" required](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
path_to_R = "/usr/local/bioinf/R/R-4.2.3/bin/"
os.environ["PATH"] = path_to_R + os.pathsep + os.environ["PATH"]

cpus = 6

# %% [markdown]
# ## Configure paths
# * adata_path: Path to anndata file
# * deseq: Path to deseq2 script
# * deseq_results: Path to results directory.
#

# %%
adata_path = "/data/projects/2023/atlas_protocol/input_data_zenodo/atlas-integrated-annotated.h5ad"
deseq = "../../bin/deseq2.R"
deseq_results = "/data/projects/2023/atlas_protocol/results/differential_expression/deseq_resdir"

# %% [markdown]
# ## Load data
#
# *anndata object

# %%
adata = sc.read_h5ad(adata_path)

# Subset adata
adata = adata[adata.obs["origin"].isin(["tumor_primary"])]
adata = adata[adata.obs["condition"].isin(["LUAD", "LUSC"])]

# %% [markdown]
# ## Pseudobulk

# %% [markdown]
# ## Get pseudobulk for entire adata

# %%
# Get pseudo-bulk profile
pdata = dc.get_pseudobulk(
    adata,
    sample_col="sample",
    groups_col="condition",
    layer="raw_counts",  # layer where raw counts are store in adata
    mode="sum",
    min_cells=0,
    min_counts=0,
)
pdata

# %% [markdown]
# ## Quality control plot
#
# From generated profile for each dataset

# %%
dc.plot_psbulk_samples(pdata, groupby=["dataset", "condition"], figsize=(13, 5))

# %% [markdown]
# Convention to filter low quality samples:
#
#     - Number of cells --> min_cells (genes  minimum total number of reads across sample)
#     - Number of counts --> min_counts (genes minimum number of counts in a number of samples)
#
# Check the frequency of genes (features) vs n. of samples and total sum of counts with  ```plot_filter_by_expr```
#

# %%
dc.plot_filter_by_expr(pdata, group="condition", min_count=10, min_total_count=15)

# %% [markdown]
# ##  Define cell type to use
#
# Specifiy for which cell type annotation level we want to run the differential expression analyses

# %%
cell_type = list(adata.obs["cell_type_coarse"].unique())

# %% [markdown]
# ## Create dictionary of adatas subsetted by cell type

# %%
adata_dict = {}
for name in cell_type:
    name_ad = name.replace(" ", "_")
    adata_name = f"{name_ad}_adata"
    adata_dict[adata_name] = adata[adata.obs["cell_type_coarse"].isin([name])]


# %% [markdown]
# ##  Function: Run DESeq2 on pseudobulk of all cell types from cell_type


# %%
def run_deseq(count_table, sample_sheet, deseq_prefix, contrast, deseq_resdir):
    """Function: Run DESeq2 on pseudobulk of all cell types from cell_type."""
    os.makedirs(deseq_resdir, exist_ok=True)

    deseq_cmd = [
        deseq,
        count_table,
        sample_sheet,
        "--cond_col",
        "condition",
        "--c1",
        contrast[0],
        "--c2",
        contrast[1],
        "--resDir",
        deseq_resdir,
        "--prefix",
        deseq_prefix,
        "--cpus",
        str(cpus),
        "--save_workspace",
    ]
    with open(deseq_resdir + "/" + deseq_prefix + ".log", "w") as stdout:
        with open(deseq_resdir + "/" + deseq_prefix + ".err", "w") as stderr:
            sp.run(deseq_cmd, capture_output=False, stdout=stdout, stderr=stderr, check=True)


# %% [markdown]
# ## Function: Homogenize all sample ids


# %%
def fix_sample_ids(pb):
    """Homogenize all sample ids."""
    repl = {}
    for k, v in dict(zip(pb.obs["condition"].index, "_" + pb.obs["condition"].values)).items():
        repl[k] = k.replace(v, "")

    return repl


# %% [markdown]
# ## Funciton: Save results from pseudobulk (samplesheet and counts) for all cell types from cell_type


# %%
def save_pseudobulk(pb, samplesheet_filename, counts_filename):
    """Save results from pseudobulk (samplesheet and counts) for all cell types from cell_type."""
    samplesheet = pb.obs.copy()
    samplesheet.reset_index(inplace=True)
    sample_ids_repl = fix_sample_ids(pb)
    bulk_df = pb.to_df().T.rename(columns=sample_ids_repl)
    bulk_df = pb.to_df().T
    bulk_df.index.name = "gene_id"
    samplesheet.to_csv(samplesheet_filename, index=False)
    bulk_df.to_csv(counts_filename)


# %% [markdown]
# ## Create pseudobulk for each celltype using the coarse cell type annotation

# %%
for ct, tmp_ad in adata_dict.items():
    pb = dc.get_pseudobulk(
        tmp_ad,
        sample_col="sample",
        groups_col="condition",
        layer="raw_counts",
        mode="sum",
        min_prop=0.05,
        min_cells=10,
        min_counts=1000,
        min_smpls=2,
    )
    if pb.obs["condition"].nunique() <= 1:
        print(f"Cell type {ct} does not have enough replicates per group")
    else:
        contrast = ["LUSC", "LUAD"]
        contrast_str = f"{contrast[0]}_vs_{contrast[1]}"
        deseq_resdir = f"{deseq_results}/{contrast_str}"

        ct = ct.replace(" ", "_")
        ct_fname = ct.replace("/", "_")
        deseq_prefix = f"{contrast_str}_{ct_fname}"

        sample_sheet = f"{deseq_results}/{deseq_prefix}.samplesheet.csv"
        count_table = f"{deseq_results}/{deseq_prefix}.counts.csv"

        save_pseudobulk(pb, sample_sheet, count_table)
        run_deseq(count_table, sample_sheet, deseq_prefix, contrast, deseq_resdir)


# %%
contrasts = [
    {"name": "LUSC_vs_LUAD", "condition": "LUSC", "reference": "LUAD"},
]
contrasts

# Cell type name without space for file name
cell_type_fn = [j.replace(" ", "_").replace("/", "_") for i, j in enumerate(cell_type)]

# %%
# Get results from folder
for contrast in contrasts:
    de_res = {}

    for ct in cell_type_fn:
        csv_file = glob.glob(deseq_resdir + "/" + contrast["name"] + "_" + ct + "_adata_DESeq2_result.tsv")
        if len(csv_file) > 0:
            res_df = pd.read_csv(csv_file[0], sep="\t")
            res_df = res_df.set_index(["gene_id"])
            ct = ct.replace("_", " ")
            # Register cell type results
            de_res[ct] = res_df
    contrast["de_res"] = de_res

# %%
# Concat and build the log2FoldChange change matrix
for contrast in contrasts:
    lfc_mat = (
        pd.concat(
            [
                res.loc[:, ["log2FoldChange"]].rename(columns={"log2FoldChange": ct})
                for ct, res in contrast["de_res"].items()
            ],
            axis=1,
            sort=True,
        )
        .fillna(0)
        .T
    )
    contrast["lfc_mat"] = lfc_mat
    display(lfc_mat)

# %%
# Concat and build the fdr
for contrast in contrasts:
    fdr_mat = (
        pd.concat(
            [res.loc[:, ["padj"]].rename(columns={"padj": ct}) for ct, res in contrast["de_res"].items()],
            axis=1,
            sort=True,
        )
        .fillna(1)
        .T
    )
    contrast["fdr_mat"] = fdr_mat
    display(fdr_mat)

# %%
# Concat and build the stat matrix
for contrast in contrasts:
    stat_mat = (
        pd.concat(
            [res.loc[:, ["stat"]].rename(columns={"stat": ct}) for ct, res in contrast["de_res"].items()],
            axis=1,
            sort=True,
        )
        .fillna(0)
        .T
    )
    contrast["stat_mat"] = stat_mat
    display(stat_mat)

# %%
logFCs = contrast["lfc_mat"]
logFCs

# %%
for contrast in contrasts:
    # Extract logFCs and pvals
    logFCs = contrast["lfc_mat"]
    pvals = contrast["fdr_mat"]

    n_sig = 0
    sig_tf = {}
    cell_type = list(logFCs.index)
    ct_dict = {}
    ct_dict["cell_type"] = []
    for ct in cell_type:
        ct_dict["cell_type"].append({"ct": ct})

    # generate a volcano plot panel for each transcription factor of interest:
    # the panels show volcano plots for each celltype in which there is a
    # signigicant transcription factor activity
    for ct in ct_dict.keys():
        n_sig = len(ct_dict[ct])

        # Calculate nrows based on ncol
        ncols = 4 if n_sig >= 4 else n_sig
        nrows = int(np.ceil(n_sig / ncols))
        # fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 4, nrows * 4))
        # empty_axs = axs.flatten()
        # axs = [{"ax": ax} for ax in zip(axs))]
        axs = ct_dict["cell_type"]

        for ax in axs:
            dc.plot_volcano(
                logFCs,
                pvals,
                ax["ct"],
                name=ax["ct"],
                top=10,
                sign_thr=0.1,
                lFCs_thr=0.5,
                return_fig=False,
                # ax=ax["ax"],
            )

        # set empty axes invisible
        # for ax in range(len(axs), len(empty_axs)):
        # empty_axs[ax].set_visible(False)

        plt.tight_layout()
        plt.show()


# %%
def plot_paired(
    adata,
    groupby,
    *,
    paired_by=None,
    var_names=None,
    show=True,
    return_fig=False,
    n_cols=4,
    panel_size=(3, 4),
    show_legend=False,
    hue=None,
    size=10,
    ylabel="expression",
    pvalues: Sequence[float] = None,
    pvalue_template=lambda x: f"unadj. p={x:.2f}, t-test",
    adjust_fdr=False,
    boxplot_properties=None,
):
    """Pairwise expression plot.Makes on panel with a paired scatterplot for each variable.

    Parameters
    ----------
    adata
        adata matrix (usually pseudobulk).
    groupby
        Column containing the grouping. Must contain exactely two different values.
    paired_by
        Column indicating the pairing (e.g. "patient")
    var_names
        Only plot these variables. Default is to plot all.
    adjust_fdr
        Adjust p-values for multiple testing using the Benjamini-Hochberg procedure.
    boxplot_properties
        Properties to pass to the boxplot function.
    hue
        Column indicating the hue.
    n_cols
        Number of columns in the figure.
    panel_size
        Size of each panel.
    pvalue_template
        Template for the p-value annotation. Must contain a single placeholder for the p-value.
    pvalues
        P-values to annotate. Must be the same length as var_names.
    return_fig
        Return the figure object.
    show
        Show the figure.
    show_legend
        Show the legend.
    size
        Size of the points.
    ylabel
        Label for the y-axis.
    """
    if boxplot_properties is None:
        boxplot_properties = {}
    groups = adata.obs[groupby].unique()
    if len(groups) != 2:
        raise ValueError("The number of groups in the group_by column must be exactely 2")

    if var_names is None:
        var_names = adata.var_names
        if len(var_names) > 20:
            warnings.warn(
                "You are plotting more than 20 variables which may be slow. "
                "Explicitly set the `var_names` paraloeter to turn this off. ",
                stacklevel=2,
            )

    X = adata[:, var_names].X
    try:
        X = X.toarray()
    except AttributeError:
        pass

    groupby_cols = [groupby]
    if paired_by is not None:
        groupby_cols.insert(0, paired_by)
    if hue is not None:
        groupby_cols.insert(0, hue)

    df = adata.obs.loc[:, groupby_cols].join(pd.DataFrame(X, index=adata.obs_names, columns=var_names))

    if paired_by is not None:
        # remove unpaired samples
        df[paired_by] = df[paired_by].astype(str)
        df.set_index(paired_by, inplace=True)
        has_matching_samples = df.groupby(paired_by).apply(lambda x: sorted(x[groupby]) == sorted(groups))
        has_matching_samples = has_matching_samples.index[has_matching_samples].values
        removed_samples = adata.obs[paired_by].nunique() - len(has_matching_samples)
        if removed_samples:
            warnings.warn(f"{removed_samples} unpaired samples removed", stacklevel=2)

        # perform statistics (paired ttest)
        if pvalues is None:
            _, pvalues = scipy.stats.ttest_rel(
                df.loc[
                    df[groupby] == groups[0],
                    var_names,
                ].loc[has_matching_samples, :],
                df.loc[
                    df[groupby] == groups[1],
                    var_names,
                ].loc[has_matching_samples],
            )

        df = df.loc[has_matching_samples, :]
        df.reset_index(drop=False, inplace=True)

    else:
        if pvalues is None:
            _, pvalues = scipy.stats.ttest_ind(
                df.loc[
                    df[groupby] == groups[0],
                    var_names,
                ],
                df.loc[
                    df[groupby] == groups[1],
                    var_names,
                ],
            )

    if adjust_fdr:
        pvalues = statsmodels.stats.multitest.fdrcorrection(pvalues)[1]

    # transform data for seaborn
    df_melt = df.melt(
        id_vars=groupby_cols,
        var_name="var",
        value_name="val",
    )

    # start plotting
    n_panels = len(var_names)
    nrows = ceil(n_panels / n_cols)
    ncols = min(n_cols, n_panels)

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(ncols * panel_size[0], nrows * panel_size[1]),
        tight_layout=True,
        squeeze=False,
    )
    axes = axes.flatten()
    if hue is None:
        hue = paired_by
    for i, (var, ax) in enumerate(zip_longest(var_names, axes)):
        if var is not None:
            sns.stripplot(
                x=groupby,
                data=df_melt.loc[lambda x: x["var"] == var],  # noqa: B023
                y="val",
                ax=ax,
                hue=hue,
                size=size,
                linewidth=1,
            )
            if paired_by is not None:
                sns.lineplot(
                    x=groupby,
                    data=df_melt.loc[lambda x: x["var"] == var],  # noqa: B023
                    hue=hue,
                    y="val",
                    ax=ax,
                    legend=False,
                    ci=None,
                )
            sns.boxplot(
                x=groupby,
                data=df_melt.loc[lambda x: x["var"] == var],  # noqa: B023
                y="val",
                ax=ax,
                color="white",
                fliersize=0,
                **boxplot_properties,
            )

            ax.set_xlabel("")
            ax.tick_params(
                axis="x",
                # rotation=0,
                labelsize=9,
            )
            ax.legend().set_visible(False)
            ax.set_ylabel(ylabel)
            ax.set_title(var + "\n" + pvalue_template(pvalues[i]))
        else:
            ax.set_visible(False)
    fig.tight_layout()

    if show_legend is True:
        axes[n_panels - 1].legend().set_visible(True)
        axes[n_panels - 1].legend(bbox_to_anchor=(1.1, 1.05))

    if show:
        plt.show()

    if return_fig:
        return fig


# %%
for ct in cell_type:
    plot_paired(
        pdata,
        "condition",
        var_names=list(contrast["de_res"][ct].index)[0:5],
        n_cols=5,
        panel_size=(2, 4),
        hue=None,
        size=10,
        ylabel="expression",
    )
