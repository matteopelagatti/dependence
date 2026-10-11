## ---------------------------------------------------------------------------
## Calibration and power of the order scan, verified from scratch.
##
## The point of the script is that nothing here is taken on trust from the
## theory.  The null table and the datasets it is applied to are independent
## draws, so the size figures are not circular; the fixed-order comparison is
## size-adjusted at its own simulated critical value, so the power comparison
## is not an artefact of one test being liberal.
##
## Runs in a few minutes.  The exactness checks -- that the profile equals
## the fixed-order statistic order by order -- are in tests/test-indepscan.R
## and run under R CMD check.
## ---------------------------------------------------------------------------

library(dependence)

N      <- 100L      # sample size
KSCAN  <- 5L        # orders scanned
RTAB   <- 20000L    # null replications per batch in the table
RNULL  <- 10000L    # fresh null datasets, for size
RPOW   <- 2000L     # datasets per alternative, for power
ALPHA  <- 0.05

WEIGHTS <- list(uniform = "uniform", geometric = "geometric")

## --------------------------------------------------------------- size -----
cat(sprintf("tabulating: n = %d, K = %d, R = %d per batch\n", N, KSCAN, RTAB))
tm <- system.time(tab <- indepscan_null(n = N, K = KSCAN, R = RTAB, seed = 101))
cat(sprintf("  %.1f s\n\n", tm[["elapsed"]]))

set.seed(999)
P <- matrix(0, RNULL, length(WEIGHTS), dimnames = list(NULL, names(WEIGHTS)))
for (i in seq_len(RNULL)) {
  x <- runif(N); y <- runif(N)
  for (j in seq_along(WEIGHTS))
    P[i, j] <- indepscan(x, y, tab, weights = WEIGHTS[[j]])$p.value
}

lev <- c(0.01, 0.05, 0.10, 0.20)
cat(sprintf("size on %d fresh null datasets (the table is an independent draw)\n",
            RNULL))
cat(sprintf("%-12s", "weights"), sprintf("%8.0f%%", 100 * lev), "\n")
for (j in seq_along(WEIGHTS))
  cat(sprintf("%-12s", names(WEIGHTS)[j]),
      sprintf("%9.4f", vapply(lev, function(a) mean(P[, j] <= a), 0)), "\n")
cat(sprintf("%-12s", "se"), sprintf("%9.4f", sqrt(lev * (1 - lev) / RNULL)), "\n")

## The p-value lives on a grid of width about 1/RTAB, so a Kolmogorov-Smirnov
## test against the CONTINUOUS uniform rejects for a reason that has nothing
## to do with calibration.  Reported as resolution, not as a verdict.
cat(sprintf("\nresolution: %d distinct p-values among %d draws (grid 1/%d)\n",
            length(unique(P[, 1])), RNULL, tab$R))

## -------------------------------------------------------------- power -----
## The fixed-order test is calibrated at its own simulated 95% point, so both
## tests sit at exactly 5% and the comparison is of power alone.
pfix <- sort(replicate(RTAB,
  indeptest(runif(N), runif(N), ties = "first")$P_stat))
cfix <- pfix[ceiling((1 - ALPHA) * RTAB)]
pord <- as.integer(max(1L, floor(N^0.31) - 1L))

alts <- list(
  monotone  = function(n) { x <- runif(n); list(x, x + 1.2 * rnorm(n)) },
  quadratic = function(n) { x <- runif(n); list(x, 4*(x-.5)^2 + .5*rnorm(n)) },
  sine      = function(n) { x <- runif(n); list(x, sin(4*pi*x) + rnorm(n)) },
  "sine 8pi" = function(n) { x <- runif(n); list(x, sin(8*pi*x) + .8*rnorm(n)) },
  circle    = function(n) { x <- runif(n); s <- sample(c(-1, 1), n, TRUE)
                            list(x, s*sqrt(pmax(0, 1-(2*x-1)^2)) + .3*rnorm(n)) },
  hetero    = function(n) { x <- abs(rnorm(n)); list(x, rnorm(n)*x + .8*rnorm(n)) }
)

cat(sprintf("\npower at %.0f%%, both tests size-adjusted; fixed order p = q = %d\n",
            100 * ALPHA, pord))
cat(sprintf("%-11s %9s %9s %9s %9s\n",
            "", "scan-unif", "scan-geom", sprintf("fixed %d", pord), "modal k"))
set.seed(4242)
for (nm in names(alts)) {
  hit <- c(0, 0, 0); km <- integer(RPOW)
  for (i in seq_len(RPOW)) {
    d <- alts[[nm]](N)
    ru <- indepscan(d[[1]], d[[2]], tab)
    rg <- indepscan(d[[1]], d[[2]], tab, weights = "geometric")
    hit <- hit + c(ru$p.value <= ALPHA, rg$p.value <= ALPHA,
                   indeptest(d[[1]], d[[2]], ties = "first")$P_stat > cfix)
    km[i] <- ru$k.min
  }
  cat(sprintf("%-11s %9.3f %9.3f %9.3f %9s\n", nm, hit[1]/RPOW, hit[2]/RPOW,
              hit[3]/RPOW, names(sort(table(km), decreasing = TRUE))[1]))
}
cat(sprintf("\nse of each entry is at most %.3f\n", 0.5 / sqrt(RPOW)))

## ------------------------------------------------------------- timing ----
cat("\none pass versus K calls to indeptest, n = 500, K = 6\n")
set.seed(7); x <- runif(500); y <- sin(4*pi*x) + .4*rnorm(500)
sp <- dependence:::scan_profile
t1 <- system.time(for (i in 1:300) sp(x, y, 6L))[["elapsed"]]
t2 <- system.time(for (i in 1:300) vapply(1:6, function(k)
  indeptest(x, y, p = k, q = k, ties = "first")$P_stat, 0))[["elapsed"]]
cat(sprintf("  scan_profile %.2f s   loop over indeptest %.2f s   %.1fx\n",
            t1, t2, t2 / t1))
