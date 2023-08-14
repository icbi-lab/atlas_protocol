# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.0
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
# In this notebook we use `decoupler` {cite}`Badia-i-Mompel2022`, a tool that contains different statistical methods to extract biological activities from omics data using prior knowledge.
# We run the method on the results from the differential gene expression analysis that compared scRNAseq data from LUAD and LUSC using celltype specific
# pseudo-bulk data.
#
# :::{note}
#
# Curated gene sets are readily available from MSigDB {cite}`Liberzon2015` and include amongst others gene
# ontology (GO)-terms {cite}`GeneOntologyConsortium2004`, KEGG {cite}`Kanehisa2000` pathways and Reactome {cite}`Gillespie2021` pathways. In addition,
# several recent efforts focus on the curation of high-quality gene signatures: PROGENy {cite}`Schubert2018`
# provides a set of reliable signatures for 14 cancer pathways derived from hundreds of perturbation
# experiments. Unlike signatures based on genes directly involved in a pathway (e.g. from KEGG),
# these signatures consist of downstream "pathway-responsive genes" that are differentially expressed
# when the pathway is perturbed (figure 1.10). This is beneficial as pathways tend to be regulated
# via post-translational modifications. Naturally, these modifications cannot be measured using
# RNA-sequencing, but requires either phosphoproteomics or targeted antibody assays. CytoSig
#  {cite}`Jiang2021` is a similar effort focussing on cytokine signaling signatures. Based on more than 2000
# cytokine treatment experiments, the authors derived 52 signatures. DoRothEA {cite}`Garcia-Alonso2019` integrates
# transcription factor activity signatures from manually curated resources, large-scale experiments
# :::
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
import os
import re
import urllib.request
from pathlib import Path

import decoupler as dc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import statsmodels.stats.multitest
from IPython.display import display

import atlas_protocol_scripts as aps

# %% [markdown]
# ## 1. Configure paths
#
# 1. Define paths to input files
#     * `adata_path`: Path to anndata file
#     * `deseq_path`: Path to directory where DESeq2 differential expression results are stored

# %%
adata_path = "../../data/input_data_zenodo/atlas-integrated-annotated.h5ad"
deseq_path = "../../data/results/differential_expression"

# %% [markdown]
# 2. Define and create output directory

# %%
results_dir = "../../data/results/10_functional_analysis"
os.makedirs(results_dir, mode=0o750, exist_ok=True)

# %% [markdown]
# 3. Define paths to network/model files. These will be automatically downloaded and stored in these files in a later step.

# %%
tfnet_file = Path(results_dir, "tf_net_dorothea_hs.tsv")
pwnet_file = Path(results_dir, "pw_net_progeny_hs.tsv")
msigdb_file = Path(results_dir, "msigdb_hs.tsv")
cytosig_file = Path(results_dir, "cytosig_signature.tsv")

# %% [markdown]
# ## 2. Load data

# %% [markdown]
# 1. Load AnnData object

# %%
adata = sc.read_h5ad(adata_path)

print(f"Anndata has: {adata.shape[0]} cells and {adata.shape[1]} genes")

# %% [markdown]
# 2. Retrieve and load DoRothEA data. We limit the analysis to the three highest confidence levels A, B and C.

# %%
if tfnet_file.exists():
    tfnet = pd.read_csv(tfnet_file, sep="\t")
else:
    tfnet = dc.get_dorothea(organism="human", levels=["A", "B", "C"])
    tfnet.to_csv(tfnet_file, sep="\t", index=False)

# %% [markdown]
# 3. Retrieve and load PROGENy data

# %%
# Retrieve Progeny db
if pwnet_file.exists():
    pwnet = pd.read_csv(pwnet_file, sep="\t")
else:
    pwnet = dc.get_progeny(organism="human", top=100)
    pwnet.to_csv(pwnet_file, sep="\t", index=False)

# %% [markdown]
# 4. Retrieve and load data from MSigDB. Here, we filter for the "hallmark" gene sets.

# %%
# Retrieve MSigDB resource
if msigdb_file.exists():
    msigdb = pd.read_csv(msigdb_file, sep="\t")
else:
    msigdb = dc.get_resource("MSigDB")

    # Filter by a desired geneset collection, for example hallmarks
    msigdb = msigdb[msigdb["collection"] == "hallmark"]
    msigdb = msigdb.drop_duplicates(["geneset", "genesymbol"])

    msigdb.to_csv(msigdb_file, sep="\t", index=False)


# %% [markdown]
# 5. Retrieve and load the CytoSig signatures

# %%
# Retrieve CytoSig signature
if cytosig_file.exists():
    cytosig_signature = pd.read_csv(cytosig_file, sep="\t")
else:
    urllib.request.urlretrieve(
        "https://github.com/data2intelligence/CytoSig/raw/master/CytoSig/signature.centroid.expand",
        cytosig_file,
    )
    cytosig_signature = pd.read_csv(cytosig_file, sep="\t")

# %% [markdown]
# ## 3. Define contrasts
#
# Here we define the gene expression contrasts/comparisons for which we to run the functional analyses.
#
# `contrasts` is a list of dicts with the following keys:
#
# * `name`: contrast name
# * `condition`: test group
# * `reference`: reference group
#
# For this tutorial, we only show one single comparison: LUAD vs. LUSC.

# %%
contrasts = [
    {"name": "LUAD_vs_LUSC", "condition": "LUAD", "reference": "LUSC"},
]
contrasts

# %% [markdown]
# 1. Create result directories for each contrast

# %%
for contrast in contrasts:
    res_dir = Path(results_dir, contrast["name"].replace(" ", "_"))
    os.makedirs(res_dir, mode=0o750, exist_ok=True)
    contrast["res_dir"] = res_dir

# %% [markdown]
# 2. Define column in `adata.obs` that contains the desired cell-type annotation. This must match the column used in the {ref}`differential_expression` section.

# %%
cell_type_class = "cell_type_coarse"

# %% [markdown]
# 3. Make list of all cell types in the specified column

# %%
cell_types = adata.obs[cell_type_class].unique()
print(f"Cell types in {cell_type_class} annotation:\n")
for ct in cell_types:
    print(ct)

# %% [markdown]
# ## 4. Read DESeq2 results

# %% [markdown]
# 1. Load TSV files created in the {ref}`differential_expression` step.

# %%
for contrast in contrasts:
    print(f"Working on: {contrast['name']}")

    de_res = {}

    for ct in cell_types:
        ct_fs = re.sub("[^0-9a-zA-Z]+", "_", ct)
        deseq_file = Path(deseq_path, contrast["name"], ct_fs, ct_fs + "_DESeq2_result.tsv")
        if os.path.exists(deseq_file):
            print(f"Reading DESeq2 result for {ct}: {deseq_file}")
            de_res[ct] = pd.read_csv(deseq_file, sep="\t").assign(cell_type=ct)
        else:
            print(f"No DESeq2 result found for: {ct}")

    contrast["cell_types"] = de_res.keys()
    contrast["de_res"] = de_res


# %% [markdown]
# 2. Convert DE results to p-value/logFC/stat matrices as used by decoupler

# %%
# Concat and build the stat matrix
for contrast in contrasts:
    lfc_mat, fdr_mat = aps.tl.long_form_df_to_decoupler(pd.concat(contrast["de_res"].values()), p_col="padj")
    stat_mat, _ = aps.tl.long_form_df_to_decoupler(
        pd.concat(contrast["de_res"].values()), log_fc_col="stat", p_col="padj"
    )
    contrast["lfc_mat"] = lfc_mat
    contrast["stat_mat"] = stat_mat
    contrast["fdr_mat"] = fdr_mat
    display(lfc_mat)
    display(fdr_mat)
    display(stat_mat)

# %% [markdown]
# ## 5. Infer pathway activities with consensus
#
# 1. Run `decoupler` consensus method to infer pathway activities from the DESeq2 result using the `PROGENy` models.\
# We use the obtained gene level `wald` statistics stored in `stat`.

# %%
# Infer pathway activities with consensus
for contrast in contrasts:
    print(contrast["name"])
    pathway_acts, pathway_pvals = dc.run_consensus(mat=contrast["stat_mat"], net=pwnet)
    contrast["pathway_acts"] = pathway_acts
    contrast["pathway_pvals"] = pathway_pvals
    display(pathway_acts)

# %% [markdown]
# 2. Generate per cell type pathway activity barplots. We use the decoupler `plot_barplot` function for plotting the pathway activities of each celltype and save the result in `png, pdf, svg` format.

# %%
# show maximum n plots in notebook, all are saved as files
show_n = 2

for contrast in contrasts:
    print(contrast["name"])

    with plt.ioff():
        p_count = 0

        for ct in contrast["cell_types"]:
            bp = dc.plot_barplot(
                contrast["pathway_acts"],
                ct,
                top=25,
                vertical=False,
                return_fig=True,
                figsize=[5, 3],
            )
            plt.title(ct)
            plt.tight_layout()

            if bp is not None:
                ct_fname = ct.replace(" ", "_").replace("/", "_")
                aps.pl.save_fig_mfmt(
                    bp,
                    res_dir=f"{contrast['res_dir']}/pathways/",
                    filename=f"{contrast['name']}_pw_acts_barplot_{ct_fname}",
                    fmt="all",
                )
                # Only show the first two plots in the notebook
                p_count += 1
                if p_count <= show_n:
                    display(bp)
                else:
                    print(f"Not showing: {ct}")

                plt.close()
            else:
                print("No plot for: " + contrast["name"] + ":" + ct)

    print(f"\nResults stored in: {contrast['res_dir']}/pathways")


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
        )
        plt.show()

        print(f"\nResults stored in: {contrast['res_dir']}/pathways")

# %% [markdown]
# ### Generate target gene expression plots for significant pathways
#
# We genereate expression plots for the target genes of pathways with significant activity differences using the DESeq2 `wald` statistics (y-axis) and the interaction `weight` (x-axis).\
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
    empty_axs = axs.flatten()
    axs = {": ".join(sig_pw.values()): {"ax": ax, "sig_pw": sig_pw} for sig_pw, ax in zip(sig_pathways, axs.flatten())}

    # Run dc.plot_targets for all significant celltype/pathway combinations using the stat values from deseq2
    for key, ax_sig_pw in axs.items():
        dc.plot_targets(
            de_res[ax_sig_pw["sig_pw"]["ct"]].set_index("gene_id"),
            stat="stat",
            source_name=ax_sig_pw["sig_pw"]["pw"],
            net=pwnet,
            top=15,
            return_fig=False,
            ax=ax_sig_pw["ax"],
        )
        ax_sig_pw["ax"].set_title(key)

    # set empty axes invisible
    for ax in range(len(axs), len(empty_axs)):
        empty_axs[ax].set_visible(False)

    plt.tight_layout()
    plt.show()

    # Save figure
    aps.pl.save_fig_mfmt(
        fig,
        res_dir=f"{contrast['res_dir']}/pathways/",
        filename=f"{contrast['name']}_pw_target_expression",
        fmt="all",
    )

    print(f"\nResults stored in: {contrast['res_dir']}/pathways")


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
# ## 6. Infer transcription factor activities with consensus
#
# Run `decoupler` consensus method to infer transcription factor activities from the DESeq2 result using the `DoRothEA` models.\
# We use the obtained gene level `wald` statistics stored in `stat`.

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

    with plt.ioff():
        p_count = 0

        for ct in contrast["cell_types"]:
            bp = dc.plot_barplot(
                contrast["tf_acts"],
                ct,
                top=25,
                vertical=False,
                return_fig=True,
                figsize=[5, 3],
            )
            plt.title(ct)
            plt.tight_layout()

            if bp is not None:
                ct_fname = ct.replace(" ", "_").replace("/", "_")
                aps.pl.save_fig_mfmt(
                    bp,
                    res_dir=f"{contrast['res_dir']}/transcription_factors/",
                    filename=f"{contrast['name']}_tf_acts_barplot_{ct_fname}",
                    fmt="all",
                )

                # Only show the first two plots in the notebook
                p_count += 1
                if p_count <= show_n:
                    display(bp)
                else:
                    print(f"Not showing: {ct}")

                plt.close()

            else:
                print("No plot for: " + contrast["name"] + ":" + ct)

        print(f"\nResults stored in: {contrast['res_dir']}/transcription_factors")


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
        )
        plt.tight_layout()
        plt.show()

    print(f"\nResults stored in: {contrast['res_dir']}/transcription_factors")

# %% [markdown]
# ### Volcano plots of expression of target genes from transcription factors of interest
#
# We genereate volcano plots for the target genes of selected transcription factors with significant activity differences using the DESeq2 `log2foldChange` (x-axis) and `padj` (y-axis) values.\
# For each transcription factor of interest a panel of celltype specific volcano plots will be created. The results are stored in `png, pdf, svg` format.
#
# `tf_of_interest`: List of transcription factors in which we are interested and for whose target genes we want to generate a volcano plot

# %%
# Define transcription factors of interest
tf_of_interest = ["CEBPA", "SOX13", "SPI1"]

for contrast in contrasts:
    # Extract logFCs and pvals
    logFCs = contrast["lfc_mat"]
    pvals = contrast["fdr_mat"]
    tf_pvals = contrast["tf_pvals"]

    # get sig ct for tfoi
    n_sig = 0
    sig_tf = {}
    for ct in contrast["cell_types"]:
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
        )

    print(f"\nResults stored in: {contrast['res_dir']}/transcription_factors")

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
# ## 7. Infer enrichment of biological terms with GSEA using significant differential expressed genes
#
# We can utilize MSigDB to assign biological terms to the differentially expressed genes. In this case, we will employ the `run_gsea` method from decoupler.

# %% [markdown]
# ### Run GSEA
# We use the wald statistics (`stat`) from the DESeq2 result to rank the genes.

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
    )
    plt.show()

    print(f"\nResults stored in: {contrast['res_dir']}/MSigDB")

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
    display(np.min(contrast["gsea_padj"].T, axis=0).sort_values().head(top_n))

# %%
for contrast in contrasts:
    gsea_norm = contrast["gsea_norm"].copy()

    # filter for enriched (positve nes score)
    gsea_norm[gsea_norm < 0] = 0

    # get top n celltypes
    top_celltypes = np.min(gsea_padj.T, axis=0).sort_values().head(top_n).index.values

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
                )
            else:
                print("No plot for: " + contrast["name"] + ":" + ct)
    plt.close()

    print(f"\nResults stored in: {contrast['res_dir']}/MSigDB")


# %% [markdown]
# ### Generate volcano plots for the most enriched term of the top n celltypes
#
# * `top_n`: number of celltypes to generate volcano plot of the most enriched term

# %%
top_n = 6

# show top n celltypes
for contrast in contrasts:
    print(contrast["name"] + "\n")
    display(np.min(contrast["gsea_padj"].T, axis=0).sort_values().tail(top_n))

# %%
for contrast in contrasts:
    logFC = contrast["lfc_mat"]
    pvals = contrast["fdr_mat"]

    # get top n celltypes
    top_celltypes = np.min(contrast["gsea_padj"].T, axis=0).sort_values().head(top_n).index.values

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
    )

    print(f"\nResults stored in: {contrast['res_dir']}/MSigDB")

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
# ## 8. CytoSig analysis
#
# We define enriched cytokine signaling signatures in the tumor cells using the *CytoSig* signature matrix an the `decoupler` consesus scoring function.

# %% [markdown]
# First we reformat the signature matrix into long format

# %%
# reformat the signature matrix
cyto_sig = pd.melt(
    cytosig_signature.rename_axis("target").reset_index(),
    var_name="source",
    id_vars=["target"],
    value_name="weight",
).reindex(columns=["source", "target", "weight"])
cyto_sig

# %% [markdown]
# Run the decoupler consensus scoring function using the DESeq2 `wald` statistics and the `CytoSig` net

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
        )
        plt.show()

        print(f"\nResults stored in: {contrast['res_dir']}/cytokine_signaling")

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
