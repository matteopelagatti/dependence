## =====================================================================
##  basis_comparison.R -- polynomial versus spline bases on the thirteen
##  schemes of the power study, for p = q = 3, 4, 5, 6, 8.
##
##  The two families are evaluated on the RANKS, as everywhere in the
##  package.  The spline basis is the package's own, bs(seq(n), df = p)
##  centred and orthonormalised: its interior knots are equally spaced on
##  the rank scale, which is the same thing as placing them at sample
##  quantiles of the data, so every knot interval holds the same number of
##  observations.  At p = 3 there are no interior knots and the spline span
##  IS the cubic polynomial span (canonical correlations 1, 1, 1), so the
##  p = 3 rows are a built-in check: the two families must agree exactly.
##
##  Size-adjusted at 5% from 2x10^4 null datasets per (n, p, basis), with a
##  size check on an independent batch.  Both B_n and P_n are recorded.
##  Every statistic is evaluated on the same datasets, so poly-vs-spline
##  differences are paired.
## =====================================================================
if (dir.exists("/home/claude/pkgtest/lib")) .libPaths(c("/home/claude/pkgtest/lib", .libPaths()))
library(dependence)
library(parallel)

SEED  <- 20261017L
NS    <- c(100L, 200L, 400L)
PS    <- c(3L, 4L, 5L, 6L, 8L)
R0    <- 20000L
R1    <- 1000L      # paired comparisons: the poly-spline differences are far more precise than a single rate
ALPHA <- 0.05
NU    <- seq(0.1, 2.0, by = 0.1)
NCORE <- if (.Platform$OS.type == "windows") 1L else max(1L, detectCores())
OUT <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  d <- if (length(m)) dirname(sub("^--file=", "", m[1L])) else getwd()
  file.path(d, "basis_results")
})
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## ------- the thirteen schemes, verbatim from 02_power_continuous.R -------
FUNCTIONAL <- c("linear", "quadratic", "cubic", "sine", "fourth root")
gen <- list(
  "linear"      = function(n, nu) { x <- runif(n); list(x = x, y = x + nu*rnorm(n)) },
  "quadratic"   = function(n, nu) { x <- runif(n); list(x = x, y = 4*(x-0.5)^2 + nu*rnorm(n)) },
  "cubic"       = function(n, nu) { x <- runif(n); z <- x - 1/3
                                    list(x = x, y = 128*z^3 - 48*z^2 - 12*z + nu*rnorm(n, 0, 3)) },
  "sine"        = function(n, nu) { x <- runif(n); list(x = x, y = sin(4*pi*x) + nu*rnorm(n, 0, 2)) },
  "fourth root" = function(n, nu) { x <- runif(n); list(x = x, y = x^0.25 + nu*rnorm(n, 0, 0.5)) },
  "circle"      = function(n, nu) { x <- runif(n); r <- rbinom(n, 1, 0.5)
                                    list(x = x, y = (2*r-1)*sqrt(pmax(0, 1-(2*x-1)^2)) + nu*rnorm(n, 0, 0.5)) },
  "two curves"  = function(n, nu) { x <- runif(n); r <- rbinom(n, 1, 0.5)
                                    list(x = x, y = 2*r*x + 0.5*(1-r)*sqrt(x) + nu*rnorm(n)) },
  "X"           = function(n, nu) { x <- runif(n); r <- rbinom(n, 1, 0.5)
                                    list(x = x, y = r*x + (1-r)*(1-x) + nu*rnorm(n, 0, 0.2)) },
  "diamond"     = function(n, nu) { x <- runif(n)
                                    y <- ifelse(x < 0.5, runif(n, 0.5-x, 0.5+x),
                                                         runif(n, x-0.5, 1.5-x))
                                    list(x = x, y = y + nu*rnorm(n, 0, 0.1)) },
  "heteroscedastic" = function(n, nu) { x <- abs(rnorm(n))
                                    list(x = x, y = rnorm(n)*x + nu*rnorm(n)) }
)
EXTRA_SCHEMES <- TRUE
if (EXTRA_SCHEMES) {
  gen[["lattice"]] <- function(n, nu) {
    ## Occupied cells: i in {1,4} pairs with j in {2,3}, i in {2,3} with
    ## j in {1,4}.  Two per row, two per column, and E[j | i] = 2.5 for all i,
    ## which is what makes every rank correlation vanish.
    i <- sample.int(4L, n, replace = TRUE)
    j <- ifelse(i == 1L | i == 4L, sample(c(2L, 3L), n, TRUE),
                                   sample(c(1L, 4L), n, TRUE))
    x <- (i - runif(n)) / 4                          # exactly uniform margins
    y <- (j - runif(n)) / 4
    list(x = x, y = y + nu*rnorm(n, 0, 0.15))
  }
  gen[["sine 8pi"]] <- function(n, nu) {
    x <- runif(n); list(x = x, y = sin(8*pi*x) + nu*rnorm(n, 0, 2))
  }
  gen[["Clayton"]] <- function(n, nu, theta = 1) {    # tau = theta/(theta+2)
    u <- runif(n); w <- runif(n)
    v <- ((w^(-theta/(1+theta)) - 1) * u^(-theta) + 1)^(-1/theta)
    list(x = u, y = v + nu*rnorm(n, 0, 0.3))
  }
  FUNCTIONAL <- c(FUNCTIONAL, "sine 8pi")            # a location relation
}
stopifnot(length(gen) == 13L)

## ------------------------------ bases --------------------------------
make_bases <- function(n) {
  B <- list(poly   = lapply(PS, function(p) dependence:::chebyshev_basis(n, p)),
            spline = lapply(PS, function(p) suppressWarnings(dependence:::ortho_spline_basis_int(n, p))))
  for (b in names(B)) for (i in seq_along(PS)) {
    M <- B[[b]][[i]]
    stopifnot(max(abs(crossprod(M)/n - diag(PS[i]))) < 1e-8, max(abs(colMeans(M))) < 1e-8)
  }
  B
}
## P_n = n ||S||_F^2 and B_n = -(n - (p+q+3)/2) sum log(1 - lambda^2), with the
## lambda the singular values of S = A'B/n, exact because A'A = B'B = n I.
all_stats <- function(rx, ry, BB, n) {
  out <- numeric(0); nm <- character(0)
  for (b in c("poly", "spline")) for (i in seq_along(PS)) {
    A <- BB[[b]][[i]]; p <- PS[i]
    S <- crossprod(A[rx, , drop = FALSE], A[ry, , drop = FALSE])/n
    l2 <- pmin(svd(S, nu = 0, nv = 0)$d^2, 1 - 1e-12)
    out <- c(out, n*sum(S^2), -(n - (2*p + 3)/2)*sum(log(1 - l2)))
    nm  <- c(nm, sprintf("P_%s_%d", b, p), sprintf("B_%s_%d", b, p))
  }
  names(out) <- nm
  out
}

null_batch <- function(n, BB, batch) {
  set.seed(SEED + 31L*n + 1000003L*batch)
  t(vapply(seq_len(R0), function(i) all_stats(sample.int(n), sample.int(n), BB, n),
           numeric(4L*length(PS))))
}
run_cell <- function(cell, BB, crit) {
  n <- cell$n
  set.seed(SEED + 1009L*n + 7919L*cell$si + 104729L*as.integer(round(cell$nu*10)))
  g <- gen[[cell$si]]
  M <- t(vapply(seq_len(R1), function(i) {
    d <- g(n, cell$nu)
    all_stats(rank(d$x, ties.method = "first"), rank(d$y, ties.method = "first"), BB, n)
  }, numeric(4L*length(PS))))
  pw <- colMeans(sweep(M, 2, crit, ">"))
  nm <- names(crit)
  data.frame(n = n, scheme = names(gen)[cell$si], nu = cell$nu,
             stat = sub("_.*", "", nm), basis = sub("^[PB]_([a-z]+)_.*", "\\1", nm),
             p = as.integer(sub(".*_", "", nm)), power = as.numeric(pw),
             stringsAsFactors = FALSE)
}

res <- list(); siz <- list()
for (n in NS) {
  BB <- make_bases(n)
  t0 <- Sys.time()
  M1 <- null_batch(n, BB, 1L); M2 <- null_batch(n, BB, 2L)
  crit <- apply(M1, 2, function(z) sort(z)[ceiling((1 - ALPHA)*R0)])
  siz[[as.character(n)]] <- data.frame(n = n, stat = names(crit),
                                       crit = as.numeric(crit),
                                       size = colMeans(sweep(M2, 2, crit, ">")),
                                       stringsAsFactors = FALSE)
  message(sprintf("n = %d: null %.0f s", n, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  cells <- list()
  for (si in seq_along(gen)) for (nu in NU) cells[[length(cells) + 1L]] <- list(n = n, si = si, nu = nu)
  t0 <- Sys.time()
  out <- mclapply(cells, run_cell, BB = BB, crit = crit, mc.cores = NCORE)
  if (any(vapply(out, function(z) !is.data.frame(z), TRUE))) stop("a cell failed")
  res[[as.character(n)]] <- do.call(rbind, out)
  message(sprintf("n = %d: %d cells %.0f s", n, length(cells),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}
write.csv(do.call(rbind, res), file.path(OUT, "basis_power.csv"), row.names = FALSE)
write.csv(do.call(rbind, siz), file.path(OUT, "basis_null.csv"),  row.names = FALSE)
message("done")
