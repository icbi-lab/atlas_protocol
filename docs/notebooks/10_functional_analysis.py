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
# # Functional analysis: progeny / dorothea

# %%
# %load_ext autoreload
# %autoreload 2

# %%
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
from IPython.display import display
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse
from threadpoolctl import threadpool_limits

import atlas_protocol_scripts as aps

warnings.simplefilter(action="ignore")
warnings.filterwarnings("ignore")

# set PATH env variable to conda env for altair_saver which is looking for npm
os.environ["PATH"] = os.path.dirname(os.environ["_"]) + os.pathsep + os.environ["PATH"]

cpus = 16
os.environ["NUMBA_NUM_THREADS"] = str(cpus)
threadpool_limits(cpus)


# %% [markdown]
# ## Configure paths

# %%
adata_path = "../../data/input_data_zenodo/atlas-integrated-annotated.h5ad"
results_dir = "../../data/results/10_functional_analysis"

tfnet_file = Path(results_dir, "tf_net_dorothea_hs.tsv")
pwnet_file = Path(results_dir, "pw_net_progeny_hs.tsv")
msigdb_file = Path(results_dir, "msigdb_hs.tsv")
cytosig_file = Path(results_dir, "cytosig_signature.tsv")


# %% [markdown]
# Create results directory

# %%
os.makedirs(results_dir, mode=0o750, exist_ok=True)

# %% [markdown]
# ## Load data

# %% [markdown]
# ### anndata object

# %%
adata = sc.read_h5ad(adata_path)

print(f"Anndata has: {adata.shape[0]} cells and {adata.shape[1]} genes")

# %% [markdown]
# ### Progeny and Dorothea data

# %%
if Path(tfnet_file).exists():
    tfnet = pd.read_csv(tfnet_file, sep="\t")
else:
    tfnet = dc.get_dorothea(organism="human", levels=["A", "B", "C"])
    tfnet.to_csv(tfnet_file, sep="\t", index=False)

# %%
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
    tfnet = pd.read_csv(msigdb_file, sep="\t")
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

# %%
contrasts = [
    {"name": "LUSC_vs_LUAD", "condition": "LUSC", "reference": "LUAD"},
]
contrasts

# %% [markdown]
# ### create result directories for each contrast

# %%
for contrast in contrasts:
    res_dir = Path(results_dir, contrast["name"].replace(" ", "_"))
    os.makedirs(res_dir, mode=0o750, exist_ok=True)
    contrast["res_dir"] = res_dir

# %% [markdown]
# ### Define cell type class to use

# %%
cell_type_class = "cell_type_coarse"

cell_types = adata.obs[cell_type_class].unique()
print(f"Cell types in {cell_type_class} annotation:")
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

# %% [markdown]
# ### Infer pathway activities with consensus

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

# %%
for contrast in contrasts:
    print(contrast["name"])
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
        else:
            print("No plot for: " + contrast["name"] + ":" + ct)


# %% [markdown]
# ### Generate pathway activity heatmap

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

# %%
# save tsv
for contrast in contrasts:
    tsv_dir = Path(contrast["res_dir"], "pathways", "tsv")
    os.makedirs(tsv_dir, mode=0o750, exist_ok=True)
    contrast["pathway_acts"].to_csv(f"{tsv_dir}/{contrast['name']}_pathway_acts.tsv", sep="\t")
    contrast["pathway_pvals"].to_csv(f"{tsv_dir}/{contrast['name']}_pathway_pvals.tsv", sep="\t")


# %% [markdown]
# ### Infer transcription factor activities with consensus

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

# %%
for contrast in contrasts:
    print(contrast["name"])
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
        else:
            print("No plot for: " + contrast["name"] + ":" + ct)


# %% [markdown]
# ### Generate transcription factor activity heatmap

# %%
# generate heatmap plot
for contrast in contrasts:
    # get acts and pvals
    pvals = contrast["tf_pvals"]
    acts = contrast["tf_acts"]

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
# ## Infer enrichment of biological terms with ORA using significant differential expressed genes

# %%
# convert deseq2 FDR results to long format
for contrast in contrasts:
    fdr_long = pd.melt(
        contrast["fdr_mat"].T.rename_axis("gene_symbol").reset_index(),
        id_vars=["gene_symbol"],
        var_name="cell_type",
        value_name="padj",
    )

    # Get top genes (FDR < 0.05)
    top_genes = fdr_long[fdr_long["padj"] < 0.05]
    contrast["top_genes"] = top_genes

# %% [markdown]
# ### Run ORA

# %%
# Infer enrichment with ora using significant deg
for contrast in contrasts:
    enr_pvals = dc.get_ora_df(
        contrast["top_genes"],
        msigdb,
        groupby="cell_type",
        features="gene_symbol",
        source="geneset",
        target="genesymbol",
    )
    contrast["enr_pvals"] = enr_pvals
    display(enr_pvals)

# %% [markdown]
# ### Genrate MSigDB heatmap from ORA scores

# %%
# Plot heatmap of ORA result
for contrast in contrasts:
    enr_pvals = contrast["enr_pvals"]
    # calculate ORA score = -log10(ora_pval)
    ora_scores = -np.log10(enr_pvals + 1e-10)

    # mark significant activities
    sig_mat = np.where(enr_pvals < 0.05, "●", "")

    # create figure
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(17, 10))

    # create heatmap
    hm = sns.heatmap(
        ora_scores,
        annot=sig_mat,
        robust=True,
        cmap="viridis",
        linewidth=0.5,
        cbar=True,
        cbar_kws={"label": "ORA score", "shrink": 0.4},
        fmt="s",
        square=True,
        xticklabels=True,
        ax=ax,
    )
    plt.tight_layout()
    aps.pl.save_fig_mfmt(
        fig,
        res_dir=f"{contrast['res_dir']}/MSigDB/",
        filename=f"{contrast['name']}_MSigDB_ORA_heatmap",
        fmt="all",
        plot_provider="mpl",
    )
    plt.show()

# %% [markdown]
# ## CytoSig analysis

# %%
cyto_sig = pd.melt(
    cytosig_signature.rename_axis("target").reset_index(), var_name="source", id_vars=["target"], value_name="weight"
).reindex(columns=["source", "target", "weight"])
cyto_sig

# %%
# Infer cytokin signaling with consensus
for contrast in contrasts:
    print(contrast["name"])

    cs_acts, cs_pvals = dc.run_consensus(mat=contrast["stat_mat"], net=cyto_sig)
    contrast["cs_acts"] = cs_acts
    contrast["cs_pvals"] = cs_pvals
    display(cs_acts)

# %%
# generate heatmap plot
for contrast in contrasts:
    # get acts and pvals
    pvals = contrast["cs_pvals"]
    acts = contrast["cs_acts"]

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

# %%
