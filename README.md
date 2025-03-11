# hambiDoubleEvo

<!-- badges: start -->
**Ecology Letters:** [![DOI:10.1111/ele.70033](https://img.shields.io/badge/DOI-10.1111/ele.70033-4F761F.svg)](https://doi.org/10.1111/ele.70033)

**Research Square:** [![DOI:10.21203/rs.3.rs-4647074/v1](https://img.shields.io/badge/DOI-10.21203/rs.3.rs--4647074/v1-B31B1B.svg)](https://doi.org/10.21203/rs.3.rs-4647074/v1)

**Code and data archive:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14017296.svg)](https://doi.org/10.5281/zenodo.14017296)
<!-- badges: end -->

[Click here to view rendered notebooks of the analysis.](https://slhogle.github.io/hambiDoubleEvo/)

## Manuscript:

‡ equal contribution, ◇ Corresponding author and equal contribution 

### Published record

**Temporal changes in the role of species sorting and evolution determine community dynamics.**\
Hoffman J<sup>‡</sup>, Hogle SL<sup>‡</sup>, Hiltunen T<sup>◇</sup>, Becks L<sup>◇</sup>. *Ecol. Lett.* (2025) [doi: 10.1111/ele.70033](https://doi.org/10.1111/ele.70033)

### Preprint

**Temporal changes in the role of species sorting and evolution determine community dynamics.**\
Hoffman J<sup>‡</sup>, Hogle SL<sup>‡</sup>, Hiltunen T<sup>◇</sup>, Becks L<sup>◇</sup>. *Research Square* (2024) [doi: 10.21203/rs.3.rs-4647074/v1](https://doi.org/10.21203/rs.3.rs-4647074/v1)

## Availability

Data and code in this GitHub repository (<https://github.com/slhogle/hambiDoubleEvo>) is provided under [GNU AGPL3](https://www.gnu.org/licenses/agpl-3.0.html). Feel free to use or remix as you see fit. 
The rendered project site is available at <https://slhogle.github.io/hambiDoubleEvo/>, which has been produced using [Quarto notebooks](https://quarto.org/). 
The content on the rendered site is released under the [CC BY 4.0.](https://creativecommons.org/licenses/by/4.0/)
This repository hosts all code and data for this project including the code necessary to fully recreate the rendered webpage.

An archived release of the code here is available from Zenodo: [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14017296.svg)](https://doi.org/10.5281/zenodo.14017296) 

Raw sequencing data using in the project is available from NCBI Bioproject [PRJNA1179357](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1179357).

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/index.html) to create reproducible environment to execute the code in this project. [See here](https://rstudio.github.io/renv/articles/renv.html#collaboration) for a brief overview on collaboration and reproduction of the entire project. To get up and running you can do:

``` r
install.packages("renv")
renv::restore()
```