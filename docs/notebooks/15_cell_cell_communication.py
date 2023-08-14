# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: atlas_protocol
#     language: python
#     name: atlas_protocol
# ---

# %% [markdown]
# # Cell-to-cell communication
#
# Effective biological functions, such as the immune response, rely heavily on the communication between cells. When cells are dissociated for single-cell sequencing, the spatial relationships between them are disrupted, leading to loss of important information. Consequently, a number of algorithms have been devised to deduce intercellular communication from gene expression data.Typically, these techniques comprise a repository of recognized ligand-receptor pairings and a computation mechanism for assessing these pairings according to their gene expression patterns.
#
# Dimitrov et al. {cite}`Dimitrov2022` conducted a comparative study of 16 ligand-receptor databases and 7 methods, which encompassed Cell-PhoneDB {cite}`Efremova2020`, CellChatDB {cite}`Jin2021`, and NATMI {cite}`Hou2020`. The LIANA package {cite}`Trei2021` was utilized to provide a unified interface to all methods, and they consolidated the ligand-receptor data from all databases into Omnipathdb. The researchers evaluated the methods and databases in relation to their correlation with spatial proximity and cytokine signaling, due to the abscence of a gold standard dataset that could be used for benchmarking purposes. Despite significant variations observed among both the methods and databases, no single technique emerged as consistently superior to the others. Methods like NicheNet {cite}`Browaeys2019`, cytotalk {cite}`Hu2021`, and SoptSC {cite}`Wang2019` go a step forward by also factoring in gene regulatory networks within a cell (intracellular communication). This enables the exploitation of the expression of downstream target genes of a receptor for predicting cell-cell interactions.
#
# Initially, the techniques outlined above were developed to evaluate intercellular communication between distinct cell types. To achieve this objective, CellPhoneDB and CellChatDB deploy a permutation test that randomizes cell type labels. While this approach is valuable in comprehending the steady state across different tissues, with the single-cell research domain now embracing perturbation experiments and extensive atlases, it has become necessary to implement methods for assessing varied patterns of intercellular communication. To takle this problem we have incorporated a straightforward yet efficient technique by considering a ligand-receptor pair as potentially perturbed if there is differential expression in at least one of the ligand and receptor. By using differentially expressed ligands/receptors (from comparisons between conditions) as input, we account for pseudoreplication bias and dropouts of lowly expressed genes.
#
# :::{seealso}
# The [Cell-cell communication](https://www.sc-best-practices.org/mechanisms/cell_cell_communication.html) chapter of the single-cell best practice book {cite}`heumosBestPracticesSinglecell2023`.
# :::

# %% [markdown]
# ## 1. Import the required libraries

# %%
import pandas as pd
import scanpy as sc

import atlas_protocol_scripts as aps

# %% [markdown]
# ## 2. Load and sanitise the input data

# %%
# Define paths
path_adata = "/home/fotakis/myScratch/atlas_tmp/atlas-integrated-annotated.h5ad"
path_ccdb = "/home/fotakis/myScratch/atlas_tmp/differential_expression/omnipathdb.tsv"
deseq2_path_prefix = "/home/fotakis/myScratch/atlas_tmp/differential_expression/"

# %%
# Load in the scRNAseq data
adata = sc.read_h5ad(path_adata)

# Load in the cellChatDB database
ccdb = pd.read_csv(path_ccdb, sep="\t")

# %%
# Keep primary sites and LUAD-LUSC conditions only
adata_primary_tumor = adata[(adata.obs["origin"] == "tumor_primary") & (adata.obs["condition"] != "NSCLC NOS")].copy()

# %%
# Sanity check
adata_primary_tumor.obs

# %% [markdown]
# ## 3. Define the Immune cell types

# %%
# Define immune cell types (to be used for ploting)
immune_cells = [
    "B cell",
    "cDC1",
    "cDC2",
    "DC mature",
    "Macrophage",
    "Macrophage alveolar",
    "Mast cell",
    "Monocyte",
    "Neutrophils",
    "NK cell",
    "pDC",
    "Plasma cell",
    "T cell CD4",
    "T cell CD8",
    "T cell regulatory",
]

# %% [markdown]
# ## 4. Perform the cell-to-cell communication analysis

# %%
# Sanity check
ccdb

# %%
# Run the analysis function (using the receptor ligand-list and the LUAD-LUSC data)
ccdba = aps.tl._cell2cell.CpdbAnalysis(
    ccdb,
    adata_primary_tumor,
    pseudobulk_group_by=["patient"],
    cell_type_column="cell_type_major",
)

# %% [markdown]
# ## 5. Use the DE results to infer cell-to-cell communications

# %%
# Load in DE results (using LUAD as reference)
de_res_tumor_cells_luad_lusc = (
    pd.read_csv(
        (deseq2_path_prefix + "/IHWallGenes.tsv").format(comparison="luad_lusc"),
        sep=",",
    )
    .fillna(1)
    .pipe(aps.tl._fdr.fdr_correction)
    .assign(group="LUAD")
)

# %%
# Sanity check
de_res_tumor_cells_luad_lusc

# %%
# Save the results
ccdb_res = ccdba.significant_interactions(de_res_tumor_cells_luad_lusc, max_pvalue=0.1)

# %% [markdown]
# ## 6. Visualisation

# %%
# Load in the results
ccdb_res = ccdba.significant_interactions(de_res_tumor_cells_luad_lusc, max_pvalue=0.1)

# %%
# Keep only the immune cells (we'll use this for plotting)
ccdb_res = ccdb_res.loc[lambda x: x["cell_type_major"].isin(immune_cells)]

# %%
# Create a list of unique top-upregulated genes
top_genes = (
    ccdb_res.loc[:, ["source_genesymbol", "fdr"]]
    .drop_duplicates()
    .sort_values("fdr")["source_genesymbol"][:30]
    .tolist()
)

# %%
# Plot the results
ccdba.plot_result(
    ccdb_res.loc[lambda x: x["source_genesymbol"].isin(top_genes)],
    title="LUAD vs LUSC: tumor cells, top 30 DE ligands",
    aggregate=False,
    cluster="heatmap",
    label_limit=80,
)
