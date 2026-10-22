## =====================================================================
##  06_timing_resampling.R
##
##  Supplement S3.6: cost of turning each statistic into a p-value.
##
##  Table S5 (results/power_continuous_scaling.csv) times one evaluation of
##  each BARE statistic at n = 1e3, 1e4, 1e5, 1e6.  This script adds:
##
##   (1) For the tests that resample each dataset (distance covariance, HHG,
##       Genest et al., MIC, BEAST): the ratio
##           time of the complete test with B = 999 replicates
##         / time of the bare statistic, as timed in Table S5,
##       measured at n = 1000.
##   (2) For the tests whose null law is simulated once per sample size and
##       then reused (order scan, Kallenberg-Ledwina, fast HHG):
##       the per-dataset time and the one-off cost of the reference law.
##   (3) B_n re-timed in the same session, so that the new rows can be
##       compared with Table S5 even if the machine load differs.
##
##  Packages: dependence, energy, HHG, MixedIndTests, minerva, BET.
##  Hoeffding's test (Hmisc::hoeffd) needs no new timing: its p-value is
##  asymptotic, so its Table S5 entries are already those of the full test.
##
##  Same data as Table S5: X standard normal, Y = X/2 + N(0,1), seed 11.
##  Run time: roughly 20-40 minutes on a laptop, mostly HHG and MIC.
##  Output: results/timing_resampling.csv
## =====================================================================

suppressPackageStartupMessages({
  library(dependence); library(energy); library(HHG)
  library(MixedIndTests); library(minerva); library(BET)
})

B_RESAMPLE <- 999L          # replicates per dataset: p-values to 0.001
N_RES      <- 1000L         # sample size for the resampling ratios
REPS_STAT  <- 3L            # medians over this many bare-statistic runs
OUT <- file.path("results", "timing_resampling.csv")

p_rule <- function(n) max(1L, as.integer(floor(n^0.31)) - 1L)
K_rule <- function(n) as.integer(floor(n^0.31)) + 1L        # two above p_rule
d_rule <- function(n) max(2L, min(4L, as.integer(floor(log2(n) / 2))))
pseudo <- function(x) rank(x, ties.method = "first") / (length(x) + 1)
dat <- function(n) { set.seed(11); x <- rnorm(n); list(x = x, y = 0.5*x + rnorm(n)) }

secs <- function(expr) {                       # elapsed seconds, NA on error
  e <- substitute(expr); env <- parent.frame()
  t0 <- proc.time()[["elapsed"]]
  ok <- tryCatch({ eval(e, env); TRUE }, error = function(err) {
    message("  error: ", conditionMessage(err)); FALSE })
  if (ok) proc.time()[["elapsed"]] - t0 else NA_real_
}
med <- function(f, reps = REPS_STAT) median(vapply(seq_len(reps), function(i) f(), 0))

rows <- list()
add <- function(test, n, what, seconds, note = "") {
  rows[[length(rows) + 1L]] <<- data.frame(test = test, n = n, what = what,
                                           secs = seconds, note = note)
  message(sprintf("%-10s n=%-8g %-28s %10.4f  %s", test, n, what, seconds, note))
}

## Defaults of the package functions, recorded for the table note.
message("BEAST defaults:");  print(formals(BET::BEAST)[c("B", "method")])
message("TestIndCopula defaults:"); print(formals(MixedIndTests::TestIndCopula)["B"])

## ---------------- (1) per-dataset resampling, n = N_RES ----------------
d <- dat(N_RES); x <- d$x; y <- d$y; n <- N_RES

## distance covariance: bare statistic as in Table S5 (dcor); permutation test
add("dCov", n, "statistic", med(function() secs(energy::dcor(x, y))))
add("dCov", n, "test", secs(energy::dcov.test(x, y, R = B_RESAMPLE)),
    paste("R =", B_RESAMPLE))

## HHG: nr.perm = 0 is the bare statistic of Table S5
hstat <- function(np) HHG::hhg.test(as.matrix(dist(x)), as.matrix(dist(y)), nr.perm = np)
add("HHG", n, "statistic", med(function() secs(hstat(0L))))
add("HHG", n, "test", secs(hstat(B_RESAMPLE)), paste("nr.perm =", B_RESAMPLE))

## Genest et al.: Table S5 used B = 2, the smallest value the package accepts
gstat <- function(B) MixedIndTests::TestIndCopula(cbind(x, y), B = B)
add("Genest", n, "statistic", med(function() secs(gstat(2L))), "B = 2")
add("Genest", n, "test", secs(gstat(B_RESAMPLE)), paste("B =", B_RESAMPLE))

## MIC: no test in minerva, so a permutation test by hand
add("MIC", n, "statistic", med(function() secs(minerva::mine_stat(x, y, measure = "mic"))))
add("MIC", n, "test", secs({
  m0 <- minerva::mine_stat(x, y, measure = "mic")
  mp <- replicate(B_RESAMPLE, minerva::mine_stat(x, sample(y), measure = "mic"))
  (1 + sum(mp >= m0)) / (B_RESAMPLE + 1) }), paste("B =", B_RESAMPLE))

## BEAST: method = "stat" is the bare statistic of Table S5; the default
## method builds the null by simulation inside the call
U <- cbind(pseudo(x), pseudo(y)); dd <- d_rule(n)
add("BEAST", n, "statistic", med(function() secs(
  BET::BEAST(U, dep = dd, unif.margin = TRUE, method = "stat"))))
add("BEAST", n, "test", secs(BET::BEAST(U, dep = dd, unif.margin = TRUE)),
    "package default")

## ------------- (2) reference law simulated once per n ------------------
## fast HHG (atoms): null table, once per n
add("HHG_fast", n, "statistic", med(function() secs(
  HHG::hhg.univariate.ind.stat(x, y, variant = "ADP-EQP-ML", nr.atoms = 40L))))
add("HHG_fast", n, "null table, once per n", secs(
  HHG::hhg.univariate.ind.nulltable(n, variant = "ADP-EQP-ML", nr.atoms = 40L)),
  "package default replicates")

## order scan and B_n at n = 1e3 ... 1e6, in this session
for (n in c(1e3, 1e4, 1e5, 1e6)) {
  n <- as.integer(n); d <- dat(n); x <- d$x; y <- d$y
  p <- p_rule(n); K <- K_rule(n)
  add("B_n", n, "statistic", med(function() secs(
    indeptest(x, y, p = p, q = p, basis = "poly", test = "Bartlett", ties = "first"))),
    paste("p = q =", p))
  add("scan", n, "profile", med(function() secs(
    dependence:::scan_profile(x, y, K, "poly", "first"))), paste("K =", K))
  if (n <= 1e4) {                              # 2e4 null profiles: minutes, not hours
    tl <- secs(law <- indepscan_null(n, K = K, R = 2e4, seed = 1))
    add("scan", n, "null law, once per n", tl, "R = 2e4")
    add("scan", n, "test given law", med(function() secs(indepscan(x, y, law))))
  }
}

res <- do.call(rbind, rows)
dir.create("results", showWarnings = FALSE)
write.csv(res, OUT, row.names = FALSE)
message("\nWritten: ", OUT)

## ratios for the table
r <- function(t) with(res, secs[test == t & what == "test"] / secs[test == t & what == "statistic"])
print(round(sapply(c("dCov", "HHG", "Genest", "MIC", "BEAST"), r), 1))
print(sessionInfo())
