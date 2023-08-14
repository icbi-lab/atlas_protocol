# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.15.0
#   kernelspec:
#     display_name: Python [conda env:conda-2023-atlas-protocol]
#     language: python
#     name: conda-env-conda-2023-atlas-protocol-py
# ---

# %% [markdown]
# # Compositional analysis

# %% [markdown]
# A central question in scRNA-seq experiments is if cell-type proportions have been changed between
# conditions. This seemingly simple question is technically challenging, due to the compositional
# nature of single-cell data. That is, if the abundance of a certain cell-type is increased, as a
# consequence, the abundance of all other cell-types will be decreased, since the overall number of
# cells profiled is limited. On top of that, cell proportions are not represented in an unbiased manner
# in scRNA-seq data, as, depending on the protocol, different cell-types are captured with different
# efficiencies {cite}`Lambrechts2018, salcherHighresolutionSinglecellAtlas2022a`.
#
# :::{note}
# Several alternative methods are available for comparing compositional data. scCODA
# {cite}`buttnerScCODABayesianModel2021` is a Bayesian model for compositional data analysis that models the uncertainty of cell-type fractions of each sample.
# It requires the definition of a reference cell-type that is assumed to
# be unchanged between conditions. tascCODA {cite}`ostnerTascCODABayesianTreeAggregated2021` is an extension of the scCODA model that additionally takes the hierarchical relationships of cell lineages into account. Propeller {cite}`phipsonPropellerTestingDifferences2022a` uses a log-linear model to model cell-type proportions and was demonstrated to have high statistical power with few biological replicates. Finally, sccomp {cite}`mangiolaRobustDifferentialComposition2022` provides a highly-flexible statistical framework that considers the presence of outliers and models group-specific variability of cell-type proportions.
#
# Another group of tools works independent of discrete cell-types and are useful for finding more
# subtle changes in functional states based on the cell × cell neighborhood graph. DA-seq {cite}`zhaoDetectionDifferentiallyAbundant2021`
# computes a differential abundance (DA)-score for each cell, based on the prevalence of conditions
# in neighborhoods of multiple sizes using a logistic regression classifier. Similarly, Milo {cite}`dannDifferentialAbundanceTesting2022` tests
# if in certain parts of the neighborhood graph cells from a certain condition are over-represented.
# Thanks to its statistics being based on a GLM with negative binomial noise model, it allows for
# flexible modeling of experimental designs and covariates.
#
# There's a [dedicated chapter](https://www.sc-best-practices.org/conditions/compositional.html) in the [single-cell best practices book](https://www.sc-best-practices.org) which provides additional information on compositional analyses.
# :::
#
# In {cite:t}`salcherHighresolutionSinglecellAtlas2022a`, we used scCODA for comparing cell-type fractions. In this section, we demonstrate how to apply it to a single-cell atlas. In this example, we are going to compare cell-type fractions between the two tumor types, *LUAD* and *LUSC*.

# %% [markdown]
# ## 1. Import the required libraries

# %%
import os

# Set tensorflow logging level to only display warnings and errors
os.environ["TF_CPP_MIN_LOG_LEVEL"] = "1"

import altair as alt
import pandas as pd
import scanpy as sc
import sccoda.util.cell_composition_data as scc_dat
import sccoda.util.comp_ana as scc_ana
import tensorflow as tf

tf.random.set_seed(0)

# %% [markdown]
# ## 2. Load the input data

# %%
adata_path = "../../data/input_data_zenodo/atlas-integrated-annotated.h5ad"

# %%
adata = sc.read_h5ad(adata_path)

# %% [markdown]
# ## 3. Compute cell-type count matrix
#
# As a first step, we compute the number of cells per sample and cell-type. This matrix is the basis for both qualitative visualization and quantitative analysis using scCODA.
#
# 1. Subset the AnnData object to the samples of interest, in our case only primary tumor samples and only patients with LUAD and LUSC.

# %%
adata_subset = adata[
    (adata.obs["origin"] == "tumor_primary") & adata.obs["condition"].isin(["LUAD", "LUSC"]),
    :,
]

# %% [markdown]
# 2. Create a DataFrame with counts per cell-type using pandas

# %%
cells_per_patient = (
    adata_subset.obs.groupby(
        # groupby needs to include all covariates of interest, the column with
        # the biological replicate (patient) and the cell-type
        ["dataset", "condition", "tumor_stage", "patient", "cell_type_coarse"],
        observed=True,
    )
    .size()
    .unstack(fill_value=0)
)

# %%
cells_per_patient

# %% [markdown]
# ## 4. Visualization as bar chart
#
# For a first (qualitative) impression, we can make a bar chart to compare cell-type fractions between conditions.
#
# 1. Transform the count matrix from above into a table of average cell-type fracitons per condition. We compute the mean cell-type fraction *per patient*, in order to give equal weight to each patient (this is particularly important if there are different numbers of cells per patient)

# %%
average_fractions_per_condition = (
    cells_per_patient.apply(lambda x: x / x.sum(), axis=1)
    .melt(ignore_index=False, value_name="frac")
    .reset_index()
    .groupby(["condition", "cell_type_coarse"], observed=True)
    .agg(mean_frac=pd.NamedAgg("frac", "mean"))
    .reset_index()
)

# %%
average_fractions_per_condition.head()

# %% [markdown]
# 2. Create a bar chart using `altair`:

# %%
alt.Chart(average_fractions_per_condition).encode(y="condition", x="mean_frac", color="cell_type_coarse").mark_bar()

# %% [markdown]
# ## 5. Quantitative analysis using scCODA
#
# 1.  Create an {class}`~anndata.AnnData` object for scCODA using {func}`~sccoda.util.cell_composition_data.from_pandas`

# %%
sccoda_data = scc_dat.from_pandas(
    cells_per_patient.reset_index(),
    # we need to specify all columns that do not contain cell-type counts as covariate columns
    covariate_columns=["patient", "dataset", "condition", "tumor_stage"],
)

# %% [markdown]
# 2. Make "condition" a categorical column. The first category is considered the "base" category (i.e. denominator of fold changes)

# %%
sccoda_data.obs["condition"] = pd.Categorical(sccoda_data.obs["condition"], categories=["LUSC", "LUAD"])

# %% [markdown]
# 3. Create the scCODA model. In `formula`, specify the condition and all covariates you would like to include into the model.
#
# :::{important}
# scCODA requires to specify a "reference cell-type" that is considered unchanged between the conditions. For studying the tumor microenvironment, we can set the reference cell-type to Epithelial cells/Tumor cells to capture changes of the tumor microenvironment
# :::

# %%
sccoda_mod = scc_ana.CompositionalAnalysis(
    sccoda_data,
    formula="condition + tumor_stage + dataset",
    # TODO consider changing reference cell-type to "tumor cells" in the final version
    reference_cell_type="Epithelial cell",
)

# %% [markdown]
# 4. Perform Hamiltonian Monte Carlo (HMC) sampling to estimate coefficients

# %%
# TODO: suitable number of iterations
sccoda_res = sccoda_mod.sample_hmc(num_results=20000)

# %% [markdown]
# 5. Set the false-discovery-rate (FDR) level for the results.
#
# :::{note}
# A smaller FDR value will produce more conservative results, but might miss some effects, while a larger FDR value selects more effects at the cost of a larger number of false discoveries. scCODA only computes a fold-change for cell-types that achieve an FDR smaller than the specified threshold. For more details, please refer to the [scCODA documentation](https://sccoda.readthedocs.io/en/latest/getting_started.html#Result-interpretation).
# :::

# %%
# TODO: suitable FDR based on final dataset
sccoda_res.set_fdr(0.5)

# %% [markdown]
# 6. Inspect the results

# %%
sccoda_res.summary()

# %% [markdown]
# 7. Compute "credible effects" based on the FDR. This will select cell-types that meet the FDR threshold.
#
# :::{note}
# `condition[T.LUAD]` refers to the coefficient of the model that captures the differences inf `LUAD` compared to `LUSC` (which se set as the reference). The name of the coefficient can be found in the summary output above. By chosing a differen coefficient (e.g. `tumor_stage[T.early]`) we could test the effects of a different variable (in this case `tumor_stage`) on the cell-type composition.
# :::

# %%
credible_effects_condition = sccoda_res.credible_effects()["condition[T.LUAD]"]

# %% [markdown]
# 8. Make a plot of log2 fold changes between LUAD and LUSC using `altair`.

# %%
alt.Chart(
    sccoda_res.effect_df.loc["condition[T.LUAD]"].loc[credible_effects_condition].reset_index(),
    title="condition",
).mark_bar().encode(
    x=alt.X("Cell Type", sort="y"),
    y="log2-fold change",
    color=alt.Color("Cell Type"),
)
