# Protocol for atlas-level integration and analysis of single cells

[![Tests][badge-tests]][link-tests]
[![Documentation][badge-docs]][link-docs]

[badge-tests]: https://img.shields.io/github/actions/workflow/status/icbi-lab/atlas_protocol/test.yaml?branch=main
[link-tests]: https://github.com/icbi-lab/atlas_protocol/actions/workflows/test.yml
[badge-docs]: https://img.shields.io/readthedocs/atlas_protocol

A computational protocol for atlas-level data integration, downstream analysis and linking with phenotypic information
from bulk RNA-seq data as performed in [Salcher et al. (2022)](<https://www.cell.com/cancer-cell/fulltext/S1535-6108(22)00499-8>).

## Getting started

This repository consists of two parts:

-   An [online version][link-docs] of the protocol
-   A python package (`atlas_protocol_scripts`) providing utility functions used in the protocol.

For **prerequisites** and instructions to **obtain data** and **software dependencies** for running
the tutorial, please see the ["before you begin"][before-you-begin] section in the documentation.

For **details on the helper functions**, see the [API documentation][link-api].

## Installation

The package with helper functions can be installed from GitHub as follows:

```bash
pip install git+https://github.com/icbi-lab/atlas_protocol.git@main
```

## Release notes

See the [changelog][changelog].

## Contact

Please use the [issue tracker][issue-tracker].

## Citation

> t.b.a

[scverse-discourse]: https://discourse.scverse.org/
[issue-tracker]: https://github.com/icbi-lab/atlas_protocol/issues
[changelog]: https://atlas_protocol.readthedocs.io/latest/changelog.html
[link-docs]: https://atlas_protocol.readthedocs.io
[link-api]: https://atlas_protocol.readthedocs.io/latest/api.html
[before-you-begin]: https://atlas-protocol.readthedocs.io/en/latest/before-you-begin.html
