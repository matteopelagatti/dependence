## =====================================================================
##  06c_timing_hhg_fast.R  --  atoms-based HHG statistic (Heller et al. 2016)
##
##  Table S5 recorded "error" at n >= 1e4.  The cause was our adapter, not
##  the method: hhg.univariate.ind.stat() defaults to
##  mmax = max(floor(sqrt(n)/2), 2), which exceeds nr.atoms = 40 once
##  n >= 6400 ("mmax should be smaller or equal to nr.atoms").
##  Here mmax = min(max(floor(sqrt(n)/2), 2), 10): the value the power study
##  used at n = 50, 100, 200 (3, 5, 7) and, from n = 400 on, the default of
##  the package's own large-sample routine, Fast.independence.test()
##  (mmax = 10, nr.atoms = 40).  Same data as Table S5.
##  Writes results/timing_hhg_fast.csv.  Run from the repo root.
## =====================================================================
dat <- function(n) { set.seed(11); x <- rnorm(n); list(x = x, y = 0.5*x + rnorm(n)) }
mmax_of <- function(n) min(max(floor(sqrt(n) / 2), 2), 10)
stat <- function(x, y) HHG::hhg.univariate.ind.stat(x, y, variant = "ADP-EQP-ML",
                         mmax = mmax_of(length(x)), nr.atoms = 40L)
el <- function(expr) { t0 <- proc.time()[["elapsed"]]; force(expr)
                       proc.time()[["elapsed"]] - t0 }
out <- NULL
for (n in as.integer(10^(3:6))) {
  d <- dat(n); x <- d$x; y <- d$y
  reps <- if (n <= 1e4) 20L else 3L
  s <- median(replicate(3, el(for (i in seq_len(reps)) stat(x, y)) / reps))
  out <- rbind(out, data.frame(n = n, mmax = mmax_of(n), what = "statistic", secs = s))
  print(tail(out, 1))
}
## null table, once per n, package defaults for the fast test (200 permutations)
for (n in as.integer(c(1e3, 1e4))) {
  s <- el(HHG::Fast.independence.test.nulltable(n, mmax = mmax_of(n), nr.atoms = 40L))
  out <- rbind(out, data.frame(n = n, mmax = mmax_of(n), what = "null table (200 perm.)", secs = s))
  print(tail(out, 1))
}
dir.create("results", showWarnings = FALSE)
write.csv(out, file.path("results", "timing_hhg_fast.csv"), row.names = FALSE)
