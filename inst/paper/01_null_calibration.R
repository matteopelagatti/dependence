## =====================================================================
##  calibrate_pq_v3.R
##
##  Largest expansion order p = q for which the chi-squared calibration of
##  the independence test keeps the empirical size close to the nominal
##  levels, for a grid of sample sizes.
##
##  Changes from v2, all of which mattered:
##
##  (a) p_max is the FIRST order that is not acceptable, minus one, not the
##      maximum over acceptable orders.  When the verdicts are not monotone
##      -- at n = 1600 v2 produced acc, acc, acc, undet, acc, undet, undet,
##      rej -- the maximum jumps over the first failure and returns 9 where
##      the intended answer is 7.
##
##  (b) Every order in the reported window is brought to a common count
##      R_FLOOR before anything else happens.  v2 escalated only where the
##      decision was in doubt, so the curve mixed orders measured at 5e4
##      with orders measured at 8.5e5; at n = 1600 the order p = 7, left at
##      the low count, fell 2.5 standard errors low on all three levels at
##      once and truncated the count.
##
##  (c) getR and getB are vectorised in p.  An environment cannot be indexed
##      by a character vector, so store[[.key(n, p_vec)]] is an error.
##
##  (d) Blocks always run to full width and the block counter is stored
##      rather than inferred from R.  v2 derived it as floor(R / R_BLOCK) and
##      permitted partial blocks, so a configuration could revisit a block
##      index, reset the RNG to a seed it had already used, and silently
##      duplicate replications -- inflating the counts and narrowing the
##      confidence intervals with no visible symptom.
##
##  (e) The cluster is stopped explicitly.  on.exit() at the top level of a
##      script registers nothing, so v2 left the workers running.
##
##  Reporting.  p_max is given at two tolerances and under two readings:
##  point, the last order whose estimated excess is within tolerance, and
##  cert, the last order whose acceptability can be demonstrated at
##  confidence CONF.  The first is what v1 measured, the second what v2
##  measured, and the gap between them is where the earlier confusion lived.
## =====================================================================

library(dependence)
library(parallel)

## ---------------------------- configuration --------------------------
N_GRID   <- c(25L, 50L, 100L, 200L, 400L, 800L, 1600L)
ALPHA    <- c(0.10, 0.05, 0.01)
TOLS     <- c(0.05, 0.10)   # tolerances reported; the largest drives the search
CONF     <- 0.95            # one-sided confidence for certification
STAT     <- "Bartlett"      # statistic driving the search; both are reported
BASIS    <- "poly"
R_FLOOR  <- 2e5             # common count every order in the window receives
R_BLOCK  <- 2e5             # width of one block; R is always a multiple of it
R_MAX    <- 1e6             # ceiling per configuration
ESC_SD   <- 3               # sharpen while the excess is within this many
                            # standard errors of the driving tolerance
WIN_DOWN <- 2L              # window around the running p_max
WIN_UP   <- 3L
P_CAP    <- 25L             # ceiling on p; a runaway guard, not a constraint
SEED     <- 20260809L
NCORES   <- max(1L, detectCores() - 1L)
OUTSTEM  <- "pq_calibration_v3"
TOL_MAIN <- max(TOLS)

## ------------------------- simulation engine -------------------------
## One worker's share of one block: nrep datasets, every order in p_vec
## evaluated on each dataset.  Returns counts [statistic, level, order].
.block <- function(nrep, n, p_vec, alpha, basis) {
  cnt <- array(0, dim = c(2L, length(alpha), length(p_vec)))
  for (i in seq_len(nrep)) {
    x <- runif(n); y <- runif(n)          # independent: H0 holds by construction
    for (j in seq_along(p_vec)) {
      tt <- indeptest(x, y, p = p_vec[j], q = p_vec[j], basis = basis,
                      test = c("Pillai", "Bartlett"))
      cnt[, , j] <- cnt[, , j] + outer(c(tt$P_pvalue, tt$B_pvalue), alpha, "<=")
    }
  }
  cnt
}

store <- new.env(parent = emptyenv())
.key  <- function(n, p) sprintf("n%d_p%d", n, p)

## Vectorised in p; each lookup is scalar because an environment cannot be
## indexed by a character vector.
getR <- function(n, p) vapply(p, function(pp) {
  z <- store[[.key(n, pp)]]; if (is.null(z)) 0 else z$R }, 0)
getB <- function(n, p) vapply(p, function(pp) {          # blocks consumed
  z <- store[[.key(n, pp)]]; if (is.null(z)) 0L else z$nblk }, 0L)

## Streams are keyed on (n, block index) and never on p, so all orders at a
## given sample size share their data and extending a configuration adds to
## its sample instead of redrawing it.  Blocks always run to full width and
## targets are rounded up to a whole number of them: a partial block would
## let two different stages collide on one stream and duplicate draws.
add_blocks <- function(cl, n, p_vec, target) {
  k <- length(cl)
  target <- ceiling(target / R_BLOCK) * R_BLOCK
  repeat {
    need <- p_vec[getR(n, p_vec) < target]
    if (!length(need)) break
    stage <- getB(n, need)
    for (s in sort(unique(stage))) {
      grp  <- need[stage == s]
      reps <- rep(R_BLOCK %/% k, k)
      if (R_BLOCK %% k)
        reps[seq_len(R_BLOCK %% k)] <- reps[seq_len(R_BLOCK %% k)] + 1L
      clusterSetRNGStream(cl, iseed = SEED + 7919 * n + 104729 * s)
      tot <- Reduce(`+`, parLapply(cl, reps, .block, n = n, p_vec = grp,
                                   alpha = ALPHA, basis = BASIS))
      for (j in seq_along(grp)) {
        kk  <- .key(n, grp[j]); cur <- store[[kk]]
        if (is.null(cur))
          cur <- list(cnt = array(0, c(2L, length(ALPHA))), R = 0, nblk = 0L)
        cur$cnt  <- cur$cnt + tot[, , j]
        cur$R    <- cur$R + sum(reps)
        cur$nblk <- cur$nblk + 1L
        store[[kk]] <- cur
      }
    }
  }
  invisible(NULL)
}

## Exact Clopper-Pearson interval, used one-sidedly at level CONF.
.ci <- function(x, R) as.numeric(
  binom.test(round(x), round(R), conf.level = 1 - 2 * (1 - CONF))$conf.int)

## Largest relative excess over the three levels, with its bounds.
excess <- function(n, p, stat = STAT) {
  cur <- store[[.key(n, p)]]
  if (is.null(cur)) stop(sprintf("no draws for n = %d, p = %d", n, p))
  i <- if (stat == "Bartlett") 2L else 1L
  e <- lo <- hi <- numeric(length(ALPHA))
  for (a in seq_along(ALPHA)) {
    b <- .ci(cur$cnt[i, a], cur$R)
    e[a]  <- cur$cnt[i, a] / cur$R / ALPHA[a] - 1
    lo[a] <- b[1] / ALPHA[a] - 1
    hi[a] <- b[2] / ALPHA[a] - 1
  }
  j <- which.max(e)
  list(e = e[j], lo = lo[j], hi = hi[j], alpha = ALPHA[j], R = cur$R,
       se = sqrt(ALPHA[j] * (1 - ALPHA[j]) / cur$R) / ALPHA[j])
}
acceptable_point <- function(n, p, tol) excess(n, p)$e  <= tol
acceptable_cert  <- function(n, p, tol) excess(n, p)$hi <= tol

## ------------------------------ search -------------------------------
cl <- makeCluster(NCORES)
invisible(clusterEvalQ(cl, library(dependence)))
clusterExport(cl, ".block")
message(sprintf("cores %d | R_floor %g | R_max %g | tolerances %s | conf %g",
                NCORES, R_FLOOR, R_MAX, paste(TOLS, collapse = ", "), CONF))

windows <- list(); p_run <- 1L
for (k in seq_along(N_GRID)) {
  n <- N_GRID[k]; cap <- min(P_CAP, n - 1L)
  lo <- max(1L, p_run - WIN_DOWN); hi <- min(cap, p_run + WIN_UP)
  add_blocks(cl, n, lo:hi, R_FLOOR)

  ## widen until the window brackets the crossing of the driving tolerance
  while (lo > 1L && !acceptable_point(n, lo, TOL_MAIN)) {
    lo <- lo - 1L; add_blocks(cl, n, lo, R_FLOOR)
  }
  while (hi < cap && acceptable_point(n, hi, TOL_MAIN)) {
    hi <- hi + 1L; add_blocks(cl, n, hi, R_FLOOR)
  }
  if (lo == 1L && !acceptable_point(n, 1L, TOL_MAIN))
    warning(sprintf("n = %d: even p = 1 exceeds the tolerance", n))
  if (hi == cap)
    warning(sprintf("n = %d: window reached the cap P_CAP = %d; raise it",
                    n, P_CAP))

  ## sharpen only the orders whose excess sits near the driving tolerance
  repeat {
    near <- Filter(function(p) {
      z <- excess(n, p)
      abs(z$e - TOL_MAIN) < ESC_SD * z$se && z$R < R_MAX
    }, lo:hi)
    if (!length(near)) break
    tgt <- min(R_MAX, max(getR(n, near)) + R_BLOCK)
    message(sprintf("    sharpening n=%d orders %s -> R=%g",
                    n, paste(near, collapse = ","), tgt))
    add_blocks(cl, n, near, tgt)
  }

  windows[[as.character(n)]] <- lo:hi
  ee <- vapply(lo:hi, function(p) excess(n, p)$e, 0)
  message(sprintf("n = %d | orders %d-%d | excess %s", n, lo, hi,
                  paste(sprintf("%+.3f", ee), collapse = " ")))
  if (is.unsorted(ee))
    message("    note: excess not monotone in p; the flat part of the curve ",
            "is not resolved at this count")

  ## carry forward the first-failure point, not the count of acceptable orders
  ok <- ee <= TOL_MAIN
  p_run <- if (all(ok)) hi else if (!ok[1]) lo else lo + which(!ok)[1] - 2L
  p_run <- max(1L, as.integer(p_run))
}

## ------------------------- first-failure rule ------------------------
pmax_of <- function(n, tol, rule = c("point", "cert")) {
  rule <- match.arg(rule); ps <- windows[[as.character(n)]]
  ok <- vapply(ps, function(p) if (rule == "point")
               acceptable_point(n, p, tol) else acceptable_cert(n, p, tol), TRUE)
  bad <- which(!ok)
  if (!length(bad))   return(max(ps))
  if (bad[1] == 1L)   return(NA_integer_)  # window did not bracket from below
  ps[bad[1] - 1L]                          # first failure, minus one
}

tab <- data.frame(n = N_GRID)
for (tol in TOLS) {
  tab[[sprintf("pmax_point_%g", tol)]] <-
    vapply(N_GRID, pmax_of, 0L, tol, "point")
  tab[[sprintf("pmax_cert_%g", tol)]]  <-
    vapply(N_GRID, pmax_of, 0L, tol, "cert")
}
print(tab, row.names = FALSE)

## ------------------------------ outputs ------------------------------
rows <- list()
for (n in N_GRID) for (p in windows[[as.character(n)]]) {
  cur <- store[[.key(n, p)]]
  for (i in 1:2) for (a in seq_along(ALPHA)) {
    b <- .ci(cur$cnt[i, a], cur$R)
    rows[[length(rows) + 1L]] <- data.frame(
      n = n, p = p, R = cur$R, stat = c("Pillai", "Bartlett")[i],
      alpha = ALPHA[a], size = cur$cnt[i, a] / cur$R, lo = b[1], hi = b[2],
      excess = cur$cnt[i, a] / cur$R / ALPHA[a] - 1,
      se = sqrt(ALPHA[a] * (1 - ALPHA[a]) / cur$R))
  }
}
curv <- do.call(rbind, rows)
curv <- curv[order(curv$n, curv$p, curv$stat, -curv$alpha), ]
write.csv(curv, paste0(OUTSTEM, "_curve.csv"), row.names = FALSE)
write.csv(tab,  paste0(OUTSTEM, ".csv"),       row.names = FALSE)

message("\nlargest excess over the three levels, ", STAT, ", by (n, p):")
for (n in N_GRID) {
  s <- vapply(windows[[as.character(n)]], function(p) {
    z <- excess(n, p); sprintf("%d:%+.3f", p, z$e) }, "")
  message(sprintf("   n=%-5d %s", n, paste(s, collapse = "  ")))
}
message("\nMonte Carlo standard error of the excess at R_floor, by level:")
for (a in ALPHA)
  message(sprintf("   alpha = %.2f : %.1f%%", a,
                  100 * sqrt(a * (1 - a) / R_FLOOR) / a))

stopCluster(cl)
message("\nWritten: ", OUTSTEM, ".csv and ", OUTSTEM, "_curve.csv")
