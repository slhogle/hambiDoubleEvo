# hambiDoubleEvo

[Click here to view rendered notebooks of the analysis.](https://slhogle.github.io/hambiDoubleEvo/)

## Manuscript:

### Published record

**Temporal changes in the role of species sorting and evolution determine community dynamics.** Hoffman J<sup>‡</sup>, Hogle SL<sup>‡</sup>, Hiltunen T<sup>◇</sup>, Becks L<sup>◇</sup>. *Ecol. Lett.* (2024/5) [doi:]()

‡ equal contribution, ◇ Corresponding author and equal contribution 

### Preprint

**Temporal changes in the role of species sorting and evolution determine community dynamics.** Hoffman J<sup>‡</sup>, Hogle SL<sup>‡</sup>, Hiltunen T<sup>◇</sup>, Becks L<sup>◇</sup>. (2024) [doi: 10.21203/rs.3.rs-4647074/v1](https://doi.org/10.1101/2023.10.31.565024)

‡ equal contribution, ◇ Corresponding author and equal contribution

## Data and Code

Data and code here is provided under GPL3. Feel free to use or remix as you see fit.

[![DOI](https://zenodo.org/badge/411568623.svg)](https://zenodo.org/badge/latestdoi/411568623)

### Project structure

1.  `/R` contains R scripts
2.  `/data` contains data that has been processed in some way for later use
3.  `/_data_raw` contains unprocessed data scraped from compute cluster
4.  `/figs` contains figures generated from R scripts

## Availability

The rendered project site is available at <https://slhogle.github.io/hambiDoubleEvo/>. The website has been produced using [Quarto notebooks](https://quarto.org/).

This GitHub repository (<https://github.com/slhogle/hambiDoubleEvo>) hosts the code and data for this project. The rendered webpage can be fully recreated using the code at <https://github.com/slhogle/hambiDoubleEvo>.

## Reproducibility

The project uses [`renv`](https://rstudio.github.io/renv/index.html) to create reproducible environment to execute the code in this project. [See here](https://rstudio.github.io/renv/articles/renv.html#collaboration) for a brief overview on collaboration and reproduction of the entire project. To get up and running you can do:

``` r
install.packages("renv")
renv::restore()
```