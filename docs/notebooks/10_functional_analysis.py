# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: 'SSH apollo-15 apollo-15: mmCD45'
#     language: ''
#     name: rik_ssh_apollo_15_apollo15mmcd45
# ---

# %% [markdown]
# # Functional analysis: PROGENy / DoRothEA / CytoSig
#
#
# To facilitate the interpretation of scRNA-seq gene expression readouts and their differences across conditions we can summarize the information
# and infer pathway activities from prior knowledge.
#
# In this notebook we use `decoupler` {cite}`Badia-i-Mompel2022` a tool that contains different statistical methods to extract biological activities from omics data using prior knowledge.
# We run the method on the results from the differential gene expression analysis that compared scRNAseq data from LUAD and LUSC using celltype specific
# pseudo-bulk data.
#
# We infer activities of the following types:
#
# * Pathways from PROGENy {cite}`Schubert2018`
# * Transcription factors from DoRothEA {cite}`Garcia-Alonso2019`
# * Molecular Signatures (MSigDB) {cite}`Liberzon2015`
# * Cytokine signaling (CytoSig) {cite}`Jiang2021`

# %%
# %load_ext autoreload
# %autoreload 2

# %%
import logging
import os
import urllib.request
import warnings
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import statsmodels.stats.multitest
from IPython.display import display
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
from threadpoolctl import threadpool_limits

import atlas_protocol_scripts as aps

# silence matplotlib logger
logging.getLogger("matplotlib").setLevel(logging.ERROR)

# silence warnings
warnings.simplefilter(action="ignore")
warnings.filterwarnings("ignore")

# set PATH env variable to conda env for altair_saver which is looking for npm
os.environ["PATH"] = os.path.dirname(os.environ["_"]) + os.pathsep + os.environ["PATH"]

cpus = 16
os.environ["NUMBA_NUM_THREADS"] = str(cpus)
threadpool_limits(cpus)


# %% [markdown]
# ## Configure paths
#
# * `adata_path`: Path to anndata file
# * `results_dir`: Path to results directory. Will be created automatically.

# %%
adata_path = "../../data/input_data_zenodo/atlas-integrated-annotated.h5ad"

results_dir = "../../data/results/10_functional_analysis"

# %% [markdown]
# Create results directory

# %%
os.makedirs(results_dir, mode=0o750, exist_ok=True)

# %%
# Path to network/model files (will be automatically generated, no need to change)
tfnet_file = Path(results_dir, "tf_net_dorothea_hs.tsv")
pwnet_file = Path(results_dir, "pw_net_progeny_hs.tsv")
msigdb_file = Path(results_dir, "msigdb_hs.tsv")
cytosig_file = Path(results_dir, "cytosig_signature.tsv")

# %% [markdown]
# ## Load data

# %% [markdown]
# ### anndata object

# %%
adata = sc.read_h5ad(adata_path)

print(f"Anndata has: {adata.shape[0]} cells and {adata.shape[1]} genes")

# %% [markdown]
# ### PROGENy and DoRothEA data

# %%
# Retrieve Dorothea db (levels A, B, C, only)
if Path(tfnet_file).exists():
    tfnet = pd.read_csv(tfnet_file, sep="\t")
else:
    tfnet = dc.get_dorothea(organism="human", levels=["A", "B", "C"])
    tfnet.to_csv(tfnet_file, sep="\t", index=False)

# %%
# Retrieve Progeny db
if Path(pwnet_file).exists():
    pwnet = pd.read_csv(pwnet_file, sep="\t")
else:
    pwnet = dc.get_progeny(organism="human", top=100)
    pwnet.to_csv(pwnet_file, sep="\t", index=False)

# %% [markdown]
# ### MSigDB data

# %%
# Retrieve MSigDB resource
if Path(msigdb_file).exists():
    msigdb = pd.read_csv(msigdb_file, sep="\t")
else:
    msigdb = dc.get_resource("MSigDB")

    # Filter by a desired geneset collection, for example hallmarks
    msigdb = msigdb[msigdb["collection"] == "hallmark"]
    msigdb = msigdb.drop_duplicates(["geneset", "genesymbol"])

    msigdb.to_csv(msigdb_file, sep="\t", index=False)


# %% [markdown]
# ### CytoSig data

# %%
# Retrieve CytoSig signature
if Path(cytosig_file).exists():
    cytosig_signature = pd.read_csv(cytosig_file, sep="\t")
else:
    urllib.request.urlretrieve(
        "https://github.com/data2intelligence/CytoSig/raw/master/CytoSig/signature.centroid.expand", cytosig_file
    )
    cytosig_signature = pd.read_csv(cytosig_file, sep="\t")

# %% [markdown]
# ## Define contrasts
#
# Here we define the gene expression contrasts/comparisons for which we to run the functional analyses.
#
# `contrasts` is a list of dicts with the following keys:
#
# * `name`: contrast name
# * `condition`: test group
# * `reference`: reference group

# %%
contrasts = [
    {"name": "LUSC_vs_LUAD", "condition": "LUSC", "reference": "LUAD"},
]
contrasts

# %% [markdown]
# ### Create result directories for each contrast

# %%
for contrast in contrasts:
    res_dir = Path(results_dir, contrast["name"].replace(" ", "_"))
    os.makedirs(res_dir, mode=0o750, exist_ok=True)
    contrast["res_dir"] = res_dir

# %% [markdown]
# ### Define cell type class to use
#
# Specifiy for which cell type annotation level we want to run the functional analyses

# %%
cell_type_class = "cell_type_coarse"

# %%
cell_types = adata.obs[cell_type_class].unique()
print(f"Cell types in {cell_type_class} annotation:\n")
for ct in cell_types:
    print(ct)

# %% [markdown]
# ### Create pseudobulk for each celltype using the coarse cell type annotation

# %%
# Store raw rounded counts in layers
adata.layers["int_counts"] = sparse.csr_matrix.ceil(adata.layers["raw_counts"])

# %%
# use decoupler to make pseudobulk
pdata = dc.get_pseudobulk(
    adata,
    sample_col="sample",
    groups_col=cell_type_class,
    layer="int_counts",
    mode="sum",
    min_cells=10,
    min_counts=1000,
)
# pdata

# %% [markdown]
# ### Run DESeq2 on pseudobulk of all celltypes from cell_type_class for each contrast

# %%
# %%capture

# Run deseq2 on pseudobulk all cell types

for contrast in contrasts:
    de_res = {}

    for ct in cell_types:
        print("Working on: " + ct)
        pb_ct = pdata[pdata.obs[cell_type_class] == ct].copy()

        # Build DESeq2 object
        dds = DeseqDataSet(
            adata=pb_ct,
            design_factors="condition",
            refit_cooks=True,
            n_cpus=cpus,
        )

        # Compute LFCs
        dds.deseq2()

        # Extract contrast between LUAD vs LUSC
        stat_res = DeseqStats(dds, contrast=["condition", contrast["condition"], contrast["reference"]], n_cpus=cpus)

        # Compute Wald test
        stat_res.summary()

        # Shrink LFCs
        coeff = "condition_" + contrast["name"]
        stat_res.lfc_shrink(coeff=coeff)

        # Register cell type results
        de_res[ct] = stat_res.results_df

    # Register results for current contrast
    contrast["de_res"] = de_res


# %% [markdown]
# Check if we got a result

# %%
contrasts[0]["de_res"]["T cell"]

# %% [markdown]
# ### Reformat and concat the deseq2 results for each contrast
#
# * Build `stat` matrix `stat_mat`

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

# %% [markdown]
# * Build `log2FoldChange` matrix `lfc_mat`

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

# %% [markdown]
# * Build `padj` matrix `fdr_mat`

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

# %% [markdown]
# ## Infer pathway activities with consensus
#
# Run `decoupler` consensus method to infer pathway activities from the DESeq2 result using the `PROGENy` models.\
# We use the obtained gene level `t-value` statistics stored in `stat`.

# %%
# Infer pathway activities with consensus
for contrast in contrasts:
    print(contrast["name"])
    pathway_acts, pathway_pvals = dc.run_consensus(mat=contrast["stat_mat"], net=pwnet)
    contrast["pathway_acts"] = pathway_acts
    contrast["pathway_pvals"] = pathway_pvals
    display(pathway_acts)

# %% [markdown]
# ### Generate per cell type pathway activity barplots
#
# We use the decoupler `plot_barplot` function for plotting the pathway activities of each celltype and save the result in `png, pdf, svg` format.

# %%
# show maximum n plots in notebook, all are saved as files
show_n = 2

for contrast in contrasts:
    print(contrast["name"])

    plt.ioff()
    p_count = 0

    for ct in cell_types:
        bp = dc.plot_barplot(contrast["pathway_acts"], ct, top=25, vertical=False, return_fig=True, figsize=[5, 3])
        plt.title(ct)
        plt.tight_layout()

        if bp is not None:
            ct_fname = ct.replace(" ", "_").replace("/", "_")
            aps.pl.save_fig_mfmt(
                bp,
                res_dir=f"{contrast['res_dir']}/pathways/",
                filename=f"{contrast['name']}_pw_acts_barplot_{ct_fname}",
                fmt="all",
                plot_provider="mpl",
            )
            # Only show the first two plots in the notebook
            p_count += 1
            if p_count <= show_n:
                display(bp)

            plt.close()
        else:
            print("No plot for: " + contrast["name"] + ":" + ct)


# %% [markdown]
# ### Generate pathway activity heatmap
#
# We use the seaborn `clustermap` function to generate a clustered heatmap of the celltype pathway activities and save the result in `png, pdf, svg` format.\
# Signigficant activity differences are marked with "●"

# %%
# generate heatmap plot
for contrast in contrasts:
    # mark significant activities
    sig_mat = np.where(contrast["pathway_pvals"] < 0.05, "●", "")

    with plt.rc_context({"figure.figsize": (5.2, 5)}):
        chm = sns.clustermap(
            contrast["pathway_acts"],
            annot=sig_mat,
            annot_kws={"fontsize": 12},
            center=0,
            cmap="bwr",
            linewidth=0.5,
            cbar_kws={"label": "Pathway activity"},
            vmin=-2,
            vmax=2,
            fmt="s",
            figsize=(7, 7),
        )
        aps.pl.reshape_clustermap(chm, cell_width=0.05, cell_height=0.05)
        aps.pl.save_fig_mfmt(
            chm,
            res_dir=f"{contrast['res_dir']}/pathways/",
            filename=f"{contrast['name']}_pw_acts_heatmap",
            fmt="all",
            plot_provider="mpl",
        )
        plt.show()

# %% [markdown]
# ### Generate target gene expression plots for significant pathways
#
# We genereate expression plots for the target genes of pathways with significant activity differences using the DESeq2 `t-values` (y-axis) and the interaction `weight` (x-axis).\
# The results are stored in `png, pdf, svg` format.

# %% [markdown]
# Get significant pathways

# %%
# filter for p-value < 0.05
for contrast in contrasts:
    p = contrast["pathway_pvals"]
    sig_pathways_idx = np.where(p < 0.05)

    # make list of celltype / pathway pairs
    sig_pathways = []

    for sig_pw in list(zip(sig_pathways_idx[0], sig_pathways_idx[1])):
        ct = p.index[sig_pw[0]]
        pw = p.columns[sig_pw[1]]
        sig_pathways.append({"ct": ct, "pw": pw})

    contrast["sig_pathways"] = sig_pathways

# %% [markdown]
# Generate plot

# %%
for contrast in contrasts:
    print(contrast["name"])
    sig_pw = contrast["sig_pathways"]

    # Calculate nrows based on ncol
    n_sig = len(sig_pw)
    ncols = 4 if n_sig >= 4 else n_sig
    nrows = int(np.ceil(n_sig / ncols))

    # Initialize the figure panel
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 4, nrows * 4))
    axs = {": ".join(sig_pw.values()): {"ax": ax, "sig_pw": sig_pw} for sig_pw, ax in zip(sig_pathways, axs.flatten())}

    # Run dc.plot_targets for all significant celltype/pathway combinations using the stat values from deseq2
    for key, ax_sig_pw in axs.items():
        dc.plot_targets(
            de_res[ax_sig_pw["sig_pw"]["ct"]],
            stat="stat",
            source_name=ax_sig_pw["sig_pw"]["pw"],
            net=pwnet,
            top=15,
            return_fig=False,
            ax=ax_sig_pw["ax"],
        )
        ax_sig_pw["ax"].set_title(key)
    plt.tight_layout()
    plt.show()

    # Save figure
    aps.pl.save_fig_mfmt(
        fig,
        res_dir=f"{contrast['res_dir']}/pathways/",
        filename=f"{contrast['name']}_pw_target_expression",
        fmt="all",
        plot_provider="mpl",
    )


# %% [markdown]
# ### Save pathway activity and p-values matrix
#
# Finally we store the pathway activity scores and p-values in `tsv` format.

# %%
# save tsv
for contrast in contrasts:
    tsv_dir = Path(contrast["res_dir"], "pathways", "tsv")
    os.makedirs(tsv_dir, mode=0o750, exist_ok=True)
    contrast["pathway_acts"].to_csv(f"{tsv_dir}/{contrast['name']}_pathway_acts.tsv", sep="\t")
    contrast["pathway_pvals"].to_csv(f"{tsv_dir}/{contrast['name']}_pathway_pvals.tsv", sep="\t")


# %% [markdown]
# ## Infer transcription factor activities with consensus
#
# Run `decoupler` consensus method to infer transcription factor activities from the DESeq2 result using the `DoRothEA` models.\
# We use the obtained gene level `t-value` statistics stored in `stat`.

# %%
# Infer transcription factor activities with consensus
for contrast in contrasts:
    print(contrast["name"])

    tf_acts, tf_pvals = dc.run_consensus(mat=contrast["stat_mat"], net=tfnet)
    contrast["tf_acts"] = tf_acts
    contrast["tf_pvals"] = tf_pvals
    display(tf_acts)

# %% [markdown]
# ### Generate per cell type transcription factor activity barplots
#
# We use the decoupler `plot_barplot` function for plotting the transcription factor activities of each celltype and save the result in `png, pdf, svg` format.

# %%
# show maximum n plots in notebook, all are saved as files
show_n = 2

for contrast in contrasts:
    print(contrast["name"])

    plt.ioff()
    p_count = 0

    for ct in cell_types:
        bp = dc.plot_barplot(contrast["tf_acts"], ct, top=25, vertical=False, return_fig=True, figsize=[5, 3])
        plt.title(ct)
        plt.tight_layout()

        if bp is not None:
            ct_fname = ct.replace(" ", "_").replace("/", "_")
            aps.pl.save_fig_mfmt(
                bp,
                res_dir=f"{contrast['res_dir']}/transcription_factors/",
                filename=f"{contrast['name']}_tf_acts_barplot_{ct_fname}",
                fmt="all",
                plot_provider="mpl",
            )

            # Only show the first two plots in the notebook
            p_count += 1
            if p_count <= show_n:
                display(bp)

            plt.close()

        else:
            print("No plot for: " + contrast["name"] + ":" + ct)


# %% [markdown]
# ### Generate transcription factor activity heatmap
#
# We use the seaborn `clustermap` function to generate a clustered heatmap of the celltype transcription factor activities and save the result in `png, pdf, svg` format.\
# Signigficant activity differences are marked with "●"

# %%
# generate heatmap plot
for contrast in contrasts:
    # get acts and pvals
    pvals = contrast["tf_pvals"].copy()
    acts = contrast["tf_acts"].copy()

    # select the columns that have significant acts (pval < 0.05)
    sig_pval_cols = pvals.columns[(pvals < 0.05).any()]

    pvals = pvals[sig_pval_cols]
    acts = acts[sig_pval_cols]

    # mark significant activities
    sig_mat = np.where(pvals < 0.05, "●", "")

    with plt.rc_context({"figure.figsize": (10, 5)}):
        chm = sns.clustermap(
            acts,
            annot=sig_mat,
            annot_kws={"fontsize": 10},
            center=0,
            cmap="bwr",
            linewidth=0.5,
            cbar_kws={"label": "TF activity"},
            vmin=-2,
            vmax=2,
            fmt="s",
            xticklabels=True,
            figsize=(4.5, 4.5),
        )
        aps.pl.reshape_clustermap(chm, cell_width=0.05, cell_height=0.05)
        aps.pl.save_fig_mfmt(
            chm,
            res_dir=f"{contrast['res_dir']}/transcription_factors/",
            filename=f"{contrast['name']}_tf_acts_heatmap",
            fmt="all",
            plot_provider="mpl",
        )
        plt.tight_layout()
        plt.show()

# %% [markdown]
# ### Volcano plots of expression of target genes from transcription factors of interest
#
# We genereate volcano plots for the target genes of selected transcription factors with significant activity differences using the DESeq2 `log2foldChange` (x-axis) and `padj` (y-axis) values.\
# For each transcription factor of interest a panel of celltype specific volcano plots will be created. The results are stored in `png, pdf, svg` format.
#
# `tf_of_interest`: List of transcription factors in which we are interested and for whose target genes we want to generate a volcano plot

# %%
# Define transcription factors of interest
tf_of_interest = ["NFKB1", "SOX2", "MYC"]

for contrast in contrasts:
    # Extract logFCs and pvals
    logFCs = contrast["lfc_mat"]
    pvals = contrast["fdr_mat"]
    tf_pvals = contrast["tf_pvals"]

    # get sig ct for tfoi
    n_sig = 0
    sig_tf = {}
    for ct in cell_types:
        for tfoi in tf_of_interest:
            if tf_pvals.loc[ct][tfoi] < 0.05:
                if tfoi not in sig_tf:
                    sig_tf[tfoi] = []
                sig_tf[tfoi].append({"ct": ct, "tf": tfoi})
                n_sig += 1

    # generate a volcano plot panel for each transcription factor of interest:
    # the panels show volcano plots for each celltype in which there is a
    # signigicant transcription factor activity
    for tf in sig_tf.keys():
        n_sig = len(sig_tf[tf])

        # Calculate nrows based on ncol
        ncols = 4 if n_sig >= 4 else n_sig
        nrows = int(np.ceil(n_sig / ncols))

        # Initialize the figure panel
        fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 4, nrows * 4))
        empty_axs = axs.flatten()
        axs = [{"ct_tf": ct_tf, "ax": ax} for ct_tf, ax in zip(sig_tf[tf], axs.flatten())]

        for ax in axs:
            dc.plot_volcano(
                logFCs,
                pvals,
                ax["ct_tf"]["ct"],
                name=ax["ct_tf"]["tf"],
                net=tfnet,
                top=10,
                sign_thr=0.1,
                lFCs_thr=0.5,
                return_fig=False,
                ax=ax["ax"],
            )

        # set empty axes invisible
        for ax in range(len(axs), len(empty_axs)):
            empty_axs[ax].set_visible(False)

        plt.tight_layout()
        plt.show()

        # Save figure
        aps.pl.save_fig_mfmt(
            fig,
            res_dir=f"{contrast['res_dir']}/transcription_factors/",
            filename=f"{contrast['name']}_{tf}_target_expression",
            fmt="all",
            plot_provider="mpl",
        )

# %% [markdown]
# ### Save transcription factor activity and p-values matrix
#
# Finally we store the transcription factor activity scores and p-values in `tsv` format.

# %%
# save tsv
for contrast in contrasts:
    tsv_dir = Path(contrast["res_dir"], "transcription_factors", "tsv")
    os.makedirs(tsv_dir, mode=0o750, exist_ok=True)
    contrast["tf_acts"].to_csv(f"{tsv_dir}/{contrast['name']}_transcription_factor_acts.tsv", sep="\t")
    contrast["tf_pvals"].to_csv(f"{tsv_dir}/{contrast['name']}_transcription_factor_pvals.tsv", sep="\t")


# %% [markdown]
# ## Infer enrichment of biological terms with GSEA using significant differential expressed genes
#
# We can utilize MSigDB to assign biological terms to the differentially expressed genes. In this case, we will employ the `run_gsea` method from decoupler.

# %% [markdown]
# ### Run GSEA
# We use the t-value (Wald test) from the DESeq2 result to rank the genes.

# %%
for contrast in contrasts:
    # run_gsea from decoupler
    gsea_estimate, gsea_norm, gsea_pvals = dc.run_gsea(
        contrast["stat_mat"],
        msigdb,
        source="geneset",
        target="genesymbol",
        times=1000,
        batch_size=10000,
        min_n=5,
        seed=4711,
        verbose=False,
    )

    contrast["gsea_estimate"] = gsea_estimate
    contrast["gsea_norm"] = gsea_norm
    contrast["gsea_pvals"] = gsea_pvals

# %% [markdown]
# Now we correct for multiple hypothesis testing using BH

# %%
for contrast in contrasts:
    # make long format
    gsea_pvals = contrast["gsea_pvals"]
    gsea_pvals_long = pd.melt(
        gsea_pvals.T.rename_axis("source").reset_index(),
        id_vars=["source"],
        var_name="cell_type",
        value_name="pval",
    )

    # run fdrcorrection and add padj to df
    gsea_pvals_long["padj"] = statsmodels.stats.multitest.fdrcorrection(gsea_pvals_long["pval"].values)[1]

    # make padj_mat
    gsea_padj = gsea_pvals_long.pivot(index="cell_type", columns="source", values="padj")
    gsea_padj.index.name = None

    # store the padj for the contrast and make sure that we preserve the order of row and columns
    contrast["gsea_padj"] = gsea_padj.reindex(index=contrast["gsea_pvals"].index).reindex(
        contrast["gsea_pvals"].columns, axis=1
    )


# %% [markdown]
# ### Generate MSigDB heatmap from GSEA normalized enrichment scores
#
# We use the seaborn `heatmap` function to generate a heatmap of the biological terms assigned by GSEA to the different celltypes and save the result in `png, pdf, svg` format.\
# Signigficantly overrepresented terms are marked with "●"
# * NES (normalized enrichment score) is color-coded

# %%
# Plot heatmap of GSEA result
for contrast in contrasts:
    gsea_padj = contrast["gsea_padj"]

    gsea_norm = contrast["gsea_norm"].copy()

    # filter for enriched (positve nes score)
    gsea_norm[gsea_norm < 0] = 0

    # mark significant enrichment
    sig_mat = np.where(gsea_padj[gsea_norm > 0] < 0.05, "●", "")

    # create figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(17, 10))

    # create heatmap
    hm = sns.heatmap(
        gsea_norm,
        annot=sig_mat,
        robust=True,
        cmap="viridis",
        linewidth=0.5,
        cbar=True,
        cbar_kws={"label": "normalized enrichment score", "shrink": 0.4},
        fmt="s",
        square=True,
        xticklabels=True,
        ax=ax,
    )
    plt.tight_layout()
    aps.pl.save_fig_mfmt(
        fig,
        res_dir=f"{contrast['res_dir']}/MSigDB/",
        filename=f"{contrast['name']}_MSigDB_GSEA_heatmap",
        fmt="all",
        plot_provider="mpl",
    )
    plt.show()

# %% [markdown]
# ### Visualize the most enriched terms as barplot
#
# We use the decoupler `plot_barplot` function to generate barplots of the most enriched biological and save the result in `png, pdf, svg` format.
# * `top_n`: top n cell types to generate barplots for. Cell types are sorted by the highest GSEA score

# %%
top_n = 4

# show top n celltypes
for contrast in contrasts:
    print(contrast["name"] + "\n")
    display(np.max(contrast["gsea_norm"].T).sort_values().tail(top_n))
    display(np.min(contrast["gsea_padj"].T).sort_values().head(top_n))

# %%
for contrast in contrasts:
    gsea_norm = contrast["gsea_norm"].copy()

    # filter for enriched (positve nes score)
    gsea_norm[gsea_norm < 0] = 0

    # get top n celltypes
    top_celltypes = np.min(gsea_padj.T).sort_values().head(top_n).index.values

    with plt.rc_context({"figure.figsize": (8, 3)}):
        for ct in top_celltypes:
            vcenter = np.max(gsea_norm.T[ct]) / 2
            bp = dc.plot_barplot(
                gsea_norm,
                ct,
                top=10,
                vertical=True,
                return_fig=True,
                vmin=0,
                vcenter=vcenter,
                cmap="Reds",
                figsize=[8, 3],
            )
            plt.title(ct)
            plt.xlabel("normalized enrichment score")
            plt.tight_layout()
            plt.show()

            if bp is not None:
                ct_fname = ct.replace(" ", "_").replace("/", "_")
                aps.pl.save_fig_mfmt(
                    bp,
                    res_dir=f"{contrast['res_dir']}/MSigDB/",
                    filename=f"{contrast['name']}_top_terms_barplot_{ct_fname}",
                    fmt="all",
                    plot_provider="mpl",
                )
            else:
                print("No plot for: " + contrast["name"] + ":" + ct)
plt.close()


# %% [markdown]
# ### Generate volcano plots for the most enriched term of the top n celltypes
#
# * `top_n`: number of celltypes to generate volcano plot of the most enriched term

# %%
top_n = 6

# %%
# show top n celltypes
np.max(contrast["gsea_norm"].T).sort_values().tail(top_n)

# %%
for contrast in contrasts:
    logFC = contrast["lfc_mat"]
    pvals = contrast["fdr_mat"]

    # get top n celltypes
    top_celltypes = np.min(contrast["gsea_padj"].T).sort_values().head(top_n).index.values

    # get top term of each celltype
    top_term = {}
    for ct in top_celltypes:
        top_term[ct] = contrast["gsea_norm"].loc[ct].idxmax()

    # Calculate nrows based on ncol
    ncols = 3 if top_n >= 3 else n_sig
    nrows = int(np.ceil(top_n / ncols))

    # Initialize the figure panel
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols * 6, nrows * 4))
    empty_axs = axs.flatten()
    axs = [[{"ct": t, "term": top_term[t]}, ax] for t, ax in zip(top_term, axs.flatten())]

    for t, ax in axs:
        dc.plot_volcano(
            logFCs,
            pvals,
            t["ct"],
            name=t["term"],
            net=msigdb,
            top=10,
            sign_thr=0.1,
            lFCs_thr=0.5,
            source="geneset",
            target="genesymbol",
            weight=None,
            return_fig=False,
            ax=ax,
        )

    # set empty axes invisible
    for ax in range(len(axs), len(empty_axs)):
        empty_axs[ax].set_visible(False)

    plt.tight_layout()
    plt.show()

    # Save figure
    aps.pl.save_fig_mfmt(
        fig,
        res_dir=f"{contrast['res_dir']}/MSigDB/",
        filename=f"{contrast['name']}_top_terms_target_expression",
        fmt="all",
        plot_provider="mpl",
    )

# %% [markdown]
# ### Save GSEA scores of enriched terms
#
# Finally we store the GSEA scores of enriched terms in `tsv` format.

# %%
# save tsv
for contrast in contrasts:
    tsv_dir = Path(contrast["res_dir"], "MSigDB", "tsv")
    os.makedirs(tsv_dir, mode=0o750, exist_ok=True)
    contrast["gsea_norm"].to_csv(f"{tsv_dir}/{contrast['name']}_MSigDB_GSEA_nes.tsv", sep="\t")
    contrast["gsea_estimate"].to_csv(f"{tsv_dir}/{contrast['name']}_MSigDB_GSEA_estimate.tsv", sep="\t")
    contrast["gsea_pvals"].to_csv(f"{tsv_dir}/{contrast['name']}_MSigDB_GSEA_pvals.tsv", sep="\t")
    contrast["gsea_padj"].to_csv(f"{tsv_dir}/{contrast['name']}_MSigDB_GSEA_padj.tsv", sep="\t")


# %% [markdown]
# ## CytoSig analysis
#
# We define enriched cytokine signaling signatures in the tumor cells using the *CytoSig* signature matrix an the `decoupler` consesus scoring function.

# %% [markdown]
# First we reformat the signature matrix into long format

# %%
# reformat the signature matrix
cyto_sig = pd.melt(
    cytosig_signature.rename_axis("target").reset_index(), var_name="source", id_vars=["target"], value_name="weight"
).reindex(columns=["source", "target", "weight"])
cyto_sig

# %% [markdown]
# Run the decoupler consensus scoring function using the DESeq2 `t-values` values and the `CytoSig` net

# %%
# Infer cytokin signaling with consensus
for contrast in contrasts:
    print(contrast["name"])

    cs_acts, cs_pvals = dc.run_consensus(mat=contrast["stat_mat"], net=cyto_sig)
    contrast["cs_acts"] = cs_acts
    contrast["cs_pvals"] = cs_pvals
    display(cs_acts)

# %% [markdown]
# ### Generate signaling activity heatmap
#
# We use the seaborn `clustermap` function to generate a clustered heatmap of the celltype signaling activities and save the result in `png, pdf, svg` format.\
# Signigficant signaling activity differences are marked with "●"

# %%
# generate heatmap plot
for contrast in contrasts:
    # get acts and pvals
    pvals = contrast["cs_pvals"].copy()
    acts = contrast["cs_acts"].copy()

    # select the columns that have significant acts (pval < 0.05)
    sig_pval_cols = pvals.columns[(pvals < 0.05).any()]

    pvals = pvals[sig_pval_cols]
    acts = acts[sig_pval_cols]

    # mark significant activities
    sig_mat = np.where(pvals < 0.05, "●", "")

    with plt.rc_context({"figure.figsize": (10, 10)}):
        chm = sns.clustermap(
            acts,
            annot=sig_mat,
            annot_kws={"fontsize": 10},
            center=0,
            cmap="bwr",
            linewidth=0.5,
            cbar_kws={"label": "signaling activity"},
            vmin=-2,
            vmax=2,
            fmt="s",
            xticklabels=True,
            figsize=(6, 6),
        )
        aps.pl.reshape_clustermap(chm, cell_width=0.05, cell_height=0.05)
        aps.pl.save_fig_mfmt(
            chm,
            res_dir=f"{contrast['res_dir']}/cytokine_signaling/",
            filename=f"{contrast['name']}_signaling_heatmap",
            fmt="all",
            plot_provider="mpl",
        )
        plt.show()

# %% [markdown]
# ### Save cytokine signaling activity scores and p-values matrix
#
# Finally we store the cytokine signaling activity scores and p-values in `tsv` format.

# %%
# save tsv
for contrast in contrasts:
    tsv_dir = Path(contrast["res_dir"], "cytokine_signaling", "tsv")
    os.makedirs(tsv_dir, mode=0o750, exist_ok=True)
    contrast["cs_acts"].to_csv(f"{tsv_dir}/{contrast['name']}_CytoSig_acts.tsv", sep="\t")
    contrast["cs_pvals"].to_csv(f"{tsv_dir}/{contrast['name']}_CytoSig_pvals.tsv", sep="\t")
