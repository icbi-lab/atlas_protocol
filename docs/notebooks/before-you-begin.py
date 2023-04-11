# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.5
#   kernelspec:
#     display_name: Python [conda env:conda-2023-atlas-protocol]
#     language: python
#     name: conda-env-conda-2023-atlas-protocol-py
# ---

# %% [markdown]
# # Before you begin

# %% [markdown]
# ## Hardware requirements
#
# For working with up to 1 million cells the absolute minimum will be a workstation with 250GB of RAM and 16 CPU cores.
# Additionally, we recommend a GPU, as data integration with scvi-tools is faster by at least one order of magnitude compared to CPU-only computing.
#
# ## Software prerequisites
# For the following instructions, we assume you are working on Linux and have set-up the [conda](https://docs.conda.io/en/latest/) package manager and a way of working with [jupyter notebooks](https://jupyter.org/).
# In case you don't, we recommend setting up [miniforge](https://github.com/conda-forge/miniforge#download) and following the jupyter lab [installation instructions](https://jupyter.org/install).
#
# ## Clone the atlas protocol repository
#
# You can obtain all notebooks and helper scripts required for this tutorial from GitHub:
#
# ```bash
# git clone https://github.com/icbi-lab/atlas_protocol.git
# cd atlas_protocol
# ```
#
# ## Installing software dependencies
#
# All required software dependencies are declared in the `environment.yml` file.
# To install all dependencies, you can create a conda environment as follows:
#
# ```bash
# conda env create -n atlas_protocol -f env/environment.yml
# conda activate atlas_protocol
# ```
#
# In order to make conda environments work with jupyter notebooks, we suggest installing [nb_conda_kernels](https://github.com/Anaconda-Platform/nb_conda_kernels).
#
# Alternatively, you can obtain a [singularity container](TODO) with all dependencies pre-installed.
#
#
# ## Obtain and preprocess single-cell datasets
#
# :::{important}
# Make sure to choose datasets carefully. Datasets may have very different characteristics, for instance
#   * generated from frozen vs. fresh tissue
#   * pre-sorted to contain only specific cell-types
#   * different sequencing platforms
#   * whole cells vs. single nuclei
#   * multimodal profiling
#
# While dataset-to-dataset differences can be mitigated by data integration, they cannot be removed completely and it is instrumental to be aware of possible biases beforehand.
# :::
#
# Publicly available single-cell datasets come in all forms and flavors. For building the atlas, we need for each dataset (1) a gene expression count matrix and (2) cell-level metadata. Gene expression data is often available from standardized repositories such as gene expression omnibus (GEO), while metadata may be available as supplementary information of the original publication. For some datasets, only read-level data can be downloaded as FASTQ files and you will need to preprocess the data from scratch. Ideally, all datasets could be re-processed from FASTQ files with a consistent reference genome and genome annotations. However, in our experience, some datasets are only available in a processed form, requiring some sort of gene identifier remapping later on.
#
# :::{seealso}
# Raw data preprocessing is beyond the scope of this tutorial. Please refer to the following sources for more details:
#
#  * The [nf-core/scrnaseq](https://nf-co.re/scrnaseq) workflow for single-cell RNA-seq preprocessing
#  * The [Raw data processing](https://www.sc-best-practices.org/introduction/raw_data_processing.html) chapter of the single-cell best practice book {cite}`heumosBestPracticesSinglecell2023`.
# :::
#
#
#
# For this tutorial, we provide four readily processed example datasets. You can download them from zenodo as follows:
#
# ```bash
# curl TODO
# ```
#
# The *Lambrechts* dataset {cite}`Lambrechts2018` has been sequenced on the 10x Genomics 3' v2 platform, the *Maynard* dataset {cite}`maynardTherapyInducedEvolutionHuman2020` using the full-length Smart-seq2 protocol {cite}`picelliSmartseq2SensitiveFulllength2013`.
# The UKIM-V dataset {cite}`salcherHighresolutionSinglecellAtlas2022a` consists of two batches, both of which have been sequenced on the *BD Rhapsody* platform. All four datasets have been generated from fresh, whole cells and have not been enriched for a particular cell-type.
#
# ## Obtain bulk RNA-seq datasets and metadata
# For the scissor analysis, bulk data needs to be prepared as an R matrix with samples in column names and gene symbols in row names containing untransformed TPM values, stored as `rds` file. The associated clinical data must be a TSV file where one column contains the sample identifiers used as rownames of the TPM matrix.
#
# For this protocol, we provide both the TPM matrix and the clinical annotation table from the TCGA LUAD and LUSC cohorts as part of the example data.  You can download them from zenodo as follows:
#
# ```bash
# curl TODO
# ```

# %%
