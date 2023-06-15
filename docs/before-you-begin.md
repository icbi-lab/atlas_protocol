# Before you begin

## Hardware requirements

For working with up to 1 million cells the absolute minimum will be a workstation with 250GB of RAM and 16 CPU cores.
Additionally, we recommend a GPU, as data integration with scvi-tools is faster by at least one order of magnitude compared to CPU-only computing.

## Software prerequisites

For the following instructions, we assume you are working on Linux and have set-up the [conda](https://docs.conda.io/en/latest/) package manager and a way of working with [jupyter notebooks](https://jupyter.org/).
In case you don't, we recommend setting up [miniforge](https://github.com/conda-forge/miniforge#download) and following the jupyter lab [installation instructions](https://jupyter.org/install).

## Clone the atlas protocol repository

You can obtain all notebooks and helper scripts required for this tutorial from GitHub:

```bash
git clone https://github.com/icbi-lab/atlas_protocol.git
cd atlas_protocol
```

## Installing software dependencies

All required software dependencies are declared in the `environment.yml` file.
To install all dependencies, you can create a conda environment as follows:

```bash
conda env create -n atlas_protocol -f env/environment.yml
conda activate atlas_protocol
# install the atlas_protocol package with helper functions
pip install git+https://github.com/icbi-lab/atlas_protocol.git
```

In order to make conda environments work with jupyter notebooks, we suggest installing [nb_conda_kernels](https://github.com/Anaconda-Platform/nb_conda_kernels).

Alternatively, you can obtain a singluarty container <!-- [singularity container](TODO) --> with all dependencies pre-installed.

## Obtain and preprocess single-cell datasets

:::{important}
Make sure to choose datasets carefully. Datasets may have very different characteristics, for instance

-   generated from frozen vs. fresh tissue
-   pre-sorted to contain only specific cell-types
-   different sequencing platforms
-   whole cells vs. single nuclei
-   multimodal profiling

While dataset-to-dataset differences can be mitigated by data integration, they cannot be removed completely and it is instrumental to be aware of possible biases beforehand.
:::

Publicly available single-cell datasets come in all forms and flavors. For building the atlas, we need for each dataset (1) a gene expression count matrix and (2) cell-level metadata. Gene expression data is often available from standardized repositories such as gene expression omnibus (GEO), while metadata may be available as supplementary information of the original publication. For some datasets, only read-level data can be downloaded as FASTQ files and you will need to preprocess the data from scratch. Ideally, all datasets could be re-processed from FASTQ files with a consistent reference genome and genome annotations. However, in our experience, some datasets are only available in a processed form, requiring some sort of gene identifier remapping later on.

:::{see also}
Raw data preprocessing is beyond the scope of this tutorial. Please refer to the following sources for more details:

-   The [nf-core/scrnaseq](https://nf-co.re/scrnaseq) workflow for single-cell RNA-seq preprocessing
-   The [Raw data processing](https://www.sc-best-practices.org/introduction/raw_data_processing.html) chapter of the single-cell best practice book {cite}`heumosBestPracticesSinglecell2023`.
-   Usefull resource to browse available datasets: [Single cell studies database](https://docs.google.com/spreadsheets/d/1En7-UV0k0laDiIfjFkdn7dggyR7jIk3WH8QgXaMOZF0/edit#gid=0)

:::

For this tutorial, we provide four readily processed example datasets. You can download them from zenodo as follows:

```bash
curl TODO
```

The _Lambrechts_ dataset {cite}`Lambrechts2018` has been sequenced on the 10x Genomics 3' v2 platform, the _Maynard_ dataset {cite}`maynardTherapyInducedEvolutionHuman2020` using the full-length Smart-seq2 protocol {cite}`picelliSmartseq2SensitiveFulllength2013`.
The UKIM-V dataset {cite}`salcherHighresolutionSinglecellAtlas2022a` consists of two batches, both of which have been sequenced on the _BD Rhapsody_ platform. All four datasets have been generated from fresh, whole cells and have not been enriched for a particular cell-type.

## Obtain bulk RNA-seq datasets and metadata

For the scissor analysis, bulk data needs to be prepared as an R matrix with samples in column names and gene symbols in row names containing untransformed TPM values, stored as `rds` file. The associated clinical data must be a TSV file where one column contains the sample identifiers used as rownames of the TPM matrix.

For this protocol, we provide both the TPM matrix and the clinical annotation table from the TCGA LUAD and LUSC cohorts as part of the example data. You can download them from zenodo as follows:

```bash
curl TODO
```

## Obtain reference genome GTF files

To facilitate integration of the four datasets, it is important to standardize the provided gene IDs. In this tutorial, we will download the GTF files from [gencode](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human)/[ensembl](https://ftp.ensembl.org/pub/) that were originally used to annotate the genes in each dataset, enabling us to remap the provided gene symbols. This remapping is necessary to resolve ambiguity in gene symbols and ensure that only counts mapped to the same genomic location are merged, using unique Ensembl IDs as identifiers.

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz
../bin/gtf_to_table.sh gencode.v32.primary_assembly.annotation.gtf.gz gencode.v32_gene_annotation_table.csv gencode
rm gencode.v32.primary_assembly.annotation.gtf.gz

wget https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
../bin/gtf_to_table.sh Homo_sapiens.GRCh38.109.gtf.gz Homo_sapiens.GRCh38.109_gene_annotation_table.csv ensembl
rm Homo_sapiens.GRCh38.109.gtf.gz
```

:::{older annotations}
`gtf_to_table.sh` will only work for newer annotations.

The first versions that do not work:
- `Gencode version 25` released 2016-07-19
- `Ensembl version 76` released 2014-07-18

For older versions use:

```bash
zcat gencode.v25.primary_assembly.annotation.gtf.gz | awk 'BEGIN{FS="\t";OFS=","}$3=="gene"{split($9,a,";");split(a[1],gene_id,"\"");split(a[4],gene_name,"\""); print gene_id[2],gene_name[2]}' | sed '1i\Geneid,GeneSymbol' > gencode.v25_gene_annotation_table.csv
```
:::