# dependence: tests of independence for variables of any kind

The **dependence** package tests the independence of two random variables,
whether they are continuous, discrete or categorical (factors in R
terminology). The test applies the same principle whatever the nature of the
two variables: it measures the canonical correlations between orthonormal
polynomial or spline expansions of the ranks (or of the level codes) of the
two variables, and summarises them with the Bartlett ($B_n$) and Pillai
($P_n$) statistics, both asymptotically chi-squared under independence.

The package provides

* `indeptest(x, y)`, the test at a given expansion order;
* `indepscan(x, y, null)`, an adaptive test that scans the expansion orders,
  with the null distribution tabulated once by `indepscan_null(n)`;
* `depanalysis(x, y)`, the deformation plots: the transformations of `x` and
  `y` that maximise their correlation, which show the form of the dependence.

## Installation

```r
# from CRAN
install.packages("dependence")

# development version from GitHub
# install.packages("remotes")
remotes::install_github("matteopelagatti/dependence", build_vignettes = TRUE)
```

## Examples

```r
library(dependence)

# numeric vs. numeric: a parabola, invisible to correlation coefficients
x <- runif(200, -1, 1)
y <- x^2 + rnorm(200, sd = 0.2)
indeptest(x, y)

# factor vs. numeric
g <- factor(sample(letters[1:3], 200, replace = TRUE))
indeptest(g, y)

# factor vs. factor: with dummy coding P_n is Pearson's chi-squared statistic
h <- factor(sample(LETTERS[1:4], 200, replace = TRUE))
indeptest(g, h, basis = "dummy")

# order scan
tab <- indepscan_null(n = 200)
indepscan(x, y, tab)

# deformation plots (require ggplot2)
m <- depanalysis(x, y)
m$plotx
m$ploty
```

See `vignette("dependence")` for a guided tour.

## Reproducibility of the paper

The scripts and results that reproduce the simulation study of the paper
are in
[`inst/paper`](https://github.com/matteopelagatti/dependence/tree/master/inst/paper)
and
[`results`](https://github.com/matteopelagatti/dependence/tree/master/results).
They are kept in
this repository only and are not part of the package distributed on CRAN.
