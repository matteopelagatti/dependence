## =====================================================================
##  06b_timing_scan.R  --  re-time the order scan after the scan_basis fix
##
##  Reinstall the package first (devtools::install() from the repo root),
##  then run from the repo root.  Appends nothing: writes
##  results/timing_scan.csv with B_n and the scan profile at
##  n = 1e3 ... 1e6, same data as Table S5.  Many repetitions at small n,
##  because the Windows timer resolves only 0.01 s.  About two minutes.
## =====================================================================
library(dependence)
p_rule <- function(n) max(1L, as.integer(floor(n^0.31)) - 1L)
K_rule <- function(n) as.integer(floor(n^0.31)) + 1L
dat <- function(n) { set.seed(11); x <- rnorm(n); list(x = x, y = 0.5*x + rnorm(n)) }
per_eval <- function(f, reps) {                 # median over 3 blocks of reps
  median(replicate(3, { t0 <- proc.time()[["elapsed"]]
                        for (i in seq_len(reps)) f()
                        (proc.time()[["elapsed"]] - t0) / reps }))
}
out <- NULL
for (j in 1:4) {
  n <- as.integer(10^(j + 2)); reps <- c(200L, 20L, 2L, 1L)[j]
  d <- dat(n); x <- d$x; y <- d$y; p <- p_rule(n); K <- K_rule(n)
  tb <- per_eval(function() indeptest(x, y, p = p, q = p, basis = "poly",
                                      test = "Bartlett", ties = "first"), reps)
  ts <- per_eval(function() dependence:::scan_profile(x, y, K, "poly", "first"), reps)
  out <- rbind(out, data.frame(n = n, p = p, K = K, B_n = tb, scan = ts))
  print(tail(out, 1))
}
dir.create("results", showWarnings = FALSE)
write.csv(out, file.path("results", "timing_scan.csv"), row.names = FALSE)
