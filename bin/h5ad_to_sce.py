#!/usr/bin/env python
"""h5ad_to_sce.py converts AnnData objects to SingleCellExperiment objects stored as rds file.

The script is based on "anndata2ri".

Usage:
    h5ad_to_sce.py INPUT_H5AD OUTPUT_SCE
"""

import anndata2ri
import rpy2.robjects as ro
import scanpy as sc
from docopt import docopt
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr


def h5ad_to_sce(input_h5ad, output_sce):
    """Convert .h5ad file to .rds file."""
    adata = sc.read_h5ad(input_h5ad)
    r_base = importr("base")
    with localconverter(anndata2ri.converter):
        sce = ro.conversion.py2rpy(adata)
    r_base.saveRDS(sce, file=output_sce, compress=False)


if __name__ == "__main__":
    arguments = docopt(__doc__)
    h5ad_to_sce(arguments["INPUT_H5AD"], arguments["OUTPUT_SCE"])
