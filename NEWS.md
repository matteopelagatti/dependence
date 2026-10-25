# dependence 1.0

* First CRAN release.
* `indeptest()`: test of independence based on the Bartlett ($B_n$) and
  Pillai ($P_n$) statistics, for numeric, integer, logical, character and
  factor variables in any combination.
* `indepscan()` and `indepscan_null()`: adaptive scan over the expansion
  order, calibrated on the simulated exact joint null distribution.
* `depanalysis()`: deformation plots of the transformations that maximise
  the correlation between the two variables.
* New vignette, `vignette("dependence")`, illustrating the three tools.
* In `indeptest()` and `depanalysis()`, the default `q = p` is now taken
  before `p` is reduced to the number of distinct values of `x` minus one,
  so each order is capped by its own variable, as documented.
