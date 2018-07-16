# topdownr

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![build status](https://travis-ci.org/sgibb/topdownr.svg?branch=master)](https://travis-ci.org/sgibb/topdownr?branch=master)
[![codecov.io](https://img.shields.io/codecov/c/github/sgibb/topdownr.svg?branch=master)](https://codecov.io/github/sgibb/topdownr/?branch=master)
[![license](http://img.shields.io/badge/license-GPL%20%28%3E=%203%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)

[![years in bioc](http://bioconductor.org/shields/years-in-bioc/topdownr.svg)](http://bioconductor.org/packages/release/bioc/html/topdownr.html)
[![bioc downloads](http://bioconductor.org/shields/downloads/topdownr.svg)](https://bioconductor.org/packages/stats/bioc/topdownr/)
Release: [![build release](http://bioconductor.org/shields/build/release/bioc/topdownr.svg)](http://bioconductor.org/checkResults/release/bioc-LATEST/topdownr/)
Devel: [![build devel](http://bioconductor.org/shields/build/devel/bioc/topdownr.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/topdownr/)

## Installation

```r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("topdownr")
```

If you want to install the development version from github
(not recommended unless you know what you are doing):

```r
if (!require("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("sgibb/topdownr")
```

## Documentation

To get started:

```r
?"topdownr-package"
vignette("data-generation", package="topdownr")
vignette("analysis", package="topdownr")
```

Development documentation: https://sgibb.github.io/topdownr/

## Questions

General questions should be asked on
the [Bioconductor support forum](https://support.bioconductor.org/),
using `topdownr` to tag the question. Feel also free to open a
GitHub [issue](https://github.com/sgibb/topdownr/issues), in
particular for bug reports.

## Contribution

See [CONTRIBUTING.md](CONTRIBUTING.md).

## Support

See [SUPPORT.md](SUPPORT.md).
