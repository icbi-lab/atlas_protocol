# Next steps

Now that you built a single-cell atlas and performed some analyses you should consider
making it available to a wider audience to maximize its utility. This section contains
our suggestions on how to share your data and code and make your analysis reproducible.

## Sharing the atlas

The [cell-x-gene](https://cellxgene.cziscience.com/) platform offers an interactive web browser
to explore single-cell datasets. You can [request](https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data)
your dataset to be added. Besides the slick web interface, cell-x-gene provides download links in `h5ad` (scverse)
and `rds` (Seurat) format. Before uploading data to cell-x-gene, the metadata must be reannotated to match standardized
ontology terms (e.g. the "cell ontology" for cell-types).

## Sharing the model

The integration with scANVI generated a pre-trained model that can be used to project additional data onto
the atlas using scArches as we have shown in scarches <!-- TODO {ref}`scarches` -->. To enable others to use this functionality
it is required that you share the scvi model and the AnnData object that was used to generate it. The pre-trained
model can be shared via the [scvi model hub](https://huggingface.co/scvi-tools) on [huggingface](https://huggingface.co/).
For more details, see the [scvi-hub upload tutorial](https://docs.scvi-tools.org/en/stable/tutorials/notebooks/scvi_hub_upload_and_large_files.html#).

## Sharing the code

Sharing the code for your entire analysis is a prerequiste for others to reproduce
your work. We recommend uploading the code on e.g. [GitHub](https://github.com).

## Sharing the environment

Single-cell analyses require an increasingly complex environment of software packages.
As results may differ slightly between different software versions, for reproducibility
it is required to declare the exact software versions used. You can export all dependencies of
in a conda environment using

```bash
conda env export > environment.yml
```

To go one step further, you can build a container (e.g. using [apptainer](https://apptainer.org/)) to obtain
a single, sharable file that also abstracts the operating system in addition to the software packages.

## A note on reproducibility

Sharing the data, code and environment is a necessary, but not sufficient condition for
reproducing an analysis {cite}`heilReproducibilityStandardsMachine2021`.

Some single cell analysis algorithms (in particular scVI/scANVI and UMAP) will yield slightly different results on
different operating systems and different hardware, trading off computational reproducibility for a significantly
faster runtime. In particular, results will differ when changing the number of cores,
or when running on a CPU/GPU of a different architecture. See also [scverse/scanpy#2014](https://github.com/scverse/scanpy/2014)
for a discussion.

To circumvent this, declare the hardware used for the analysis and consider sharing intermediate results,
such as the scVI and UMAP embeddings.
