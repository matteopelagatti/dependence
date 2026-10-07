## =====================================================================
##  02_power_continuous.R   (resumable)
##
##  Size-adjusted power of B_n against its competitors on the ten continuous
##  dependence schemes, at n = 50, 100, 200.
##
##  WHY THIS VERSION.  The first draft built 900 cells and calibrated all
##  seven statistics separately in every one of them, with the n = 200 cells
##  -- roughly three quarters of the total work, the leading statistics
##  being quadratic in n -- queued last, and no way to see progress or to
##  resume after an interruption.  Three changes:
##
##  (a) RESUMABLE.  Each cell is written to results/ as soon as it finishes
##      and is skipped on a re-run.  Cells are self-seeded, so this changes
##      nothing about the answers.  Kill and restart freely.
##
##  (b) RANK-INVARIANT STATISTICS ARE CALIBRATED ONCE PER (n, p, q), not
##      once per cell.  A statistic that is a function of the ranks alone is
##      distribution-free under independence, so its null distribution does
##      not depend on the scheme or the noise level.  Which statistics these
##      are is DETECTED at startup rather than assumed, by checking whether
##      the statistic is unchanged when each variable is replaced by its
##      ranks.  Expect B_n, P_n, Hoeffding and the Genest statistic, which is
##      a functional of the empirical copula, to qualify; distance covariance
##      and HHG use the values and will not.
##
##  (c) EXPENSIVE CELLS FIRST, which is the right order for load-balanced
##      scheduling and means the costly n = 200 results arrive early rather
##      than last.
##
##  Together these cut the work by roughly forty per cent, on top of the
##  thousandfold already saved by computing statistics rather than calling
##  each package's own resampling test.
##
## ---------------------------------------------------------------------
##  REVISION, October 2026.  The four comparisons asked for by the JRSS-B
##  associate editor: BET and BEAST (Zhang and co-authors), Chatterjee's xi,
##  and the atoms-based fast variant of HHG.  Four supporting changes.
##
##  (d) CACHE INVALIDATION, and read this before re-running.  Cells and
##      critical-value tables from the previous study were written with the
##      OLD set of statistics.  Resuming naively would skip those cells and
##      silently leave the four new statistics missing from most of the
##      grid.  A cached file is therefore accepted only if it contains every
##      statistic currently being kept; otherwise it is recomputed.  Expect
##      a full re-run the first time, since no old cell has BET in it.
##
##  (e) PER-STATISTIC EVALUATION CAPS.  BEAST resamples internally to build
##      its adaptive weights, so one evaluation can cost orders of magnitude
##      more than any other statistic here.  EVAL_CAP lets a named statistic
##      use fewer replications than the rest; it is evaluated on a PREFIX of
##      the same datasets, so all tests still see identical data and only
##      its standard error grows, which the se column reflects.  Set the
##      caps from the timing table printed at startup, not from guesswork.
##
##  (f) COST PROJECTION BEFORE COMMITMENT.  Every adapter is timed and the
##      projected single-core hours printed per statistic, with a hard stop
##      above MAX_PROJECTED_HOURS.  A study that would take three weeks
##      should be discovered in the first minute.
##
##  (g) TOLERANT ACCESSORS.  The new packages return test objects whose
##      field names have changed across versions, and this harness wants a
##      statistic.  Each new adapter tries a list of candidate fields
##      through pluck1(); the self-test turns a wrong guess into a visible
##      failure at startup rather than a corrupted cell.  Verify each
##      against the package's own test function on one dataset before the
##      full run, as was done for dCov, HHG and Genest.
## =====================================================================

library(dependence)
library(parallel)

## ---------------------------- configuration --------------------------
N_GRID  <- c(50L, 100L, 200L)
NU      <- seq(0.1, 2.0, by = 0.1)
ALPHA   <- 0.05
R1      <- 1e4        # alternative datasets per cell
R0      <- 5e3        # null datasets per cell, for statistics that need them
R0_DF   <- 2e4        # null datasets per (n, p, q) for rank-invariant ones
USE_MIC <- TRUE
USE_NEW <- TRUE       # BET, BEAST, xi, fast HHG: the associate editor's four
SECONDARY_Q1 <- FALSE # run the primary study first; set TRUE for the q = 1
                      # variant on the five functional schemes afterwards
SEED    <- 20260810L
NCORES  <- max(1L, detectCores() - 1L)
OUTDIR  <- "results"
OUTSTEM <- "power_continuous"
dir.create(OUTDIR, showWarnings = FALSE)

## Replications for named statistics where the default is too expensive.
## Empty means everything runs at full length.  Fill from the timing table.
## Example, if BEAST costs 40 ms an evaluation:  EVAL_CAP <- list(BEAST = 2e3)
EVAL_CAP <- list()
MAX_PROJECTED_HOURS <- 400   # hard stop, single-core hours

p_rule <- function(n) max(1L, as.integer(floor(n^0.31)) - 1L)
FUNCTIONAL <- c("linear", "quadratic", "cubic", "sine", "fourth root")

## Binary expansion depth for BET and BEAST: the analogue of rule (12).
## Resolution should grow with n but stay coarse enough that the 2^d x 2^d
## cells are not mostly empty.  Gives 2, 3, 3 at n = 50, 100, 200, against
## p = q = 2, 3, 4 for our own test, so the two look at a comparable number
## of directions and neither is favoured by the tuning.
d_rule <- function(n) max(2L, min(4L, as.integer(floor(log2(n) / 2))))

cap_for <- function(nm, default)
  if (is.null(EVAL_CAP[[nm]])) default else min(default, as.integer(EVAL_CAP[[nm]]))

## Some statistics consume random numbers of their own: BEAST resamples to
## build its adaptive weights, and the Genest adapter runs a bootstrap at
## B = 2.  Evaluated inside the replication loop they advance the generator,
## so the sequence of DATASETS a cell sees would depend on which statistics
## are being kept -- and therefore change if an adapter is added, dropped by
## the self-test, or capped by EVAL_CAP.  That silently breaks both the
## identical-data property and reproducibility across runs.  Freezing the
## generator around each evaluation block makes the data depend only on the
## cell's own seed.  This also means a future run that adds a statistic can
## regenerate exactly the same datasets, which is what would make an
## incremental extension of an existing results/ directory possible; it is
## not possible for the cells already on disk, which were written before
## this fix and whose data stream is entangled with the Genest bootstrap.
with_frozen_rng <- function(expr) {
  st <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(assign(".Random.seed", st, envir = .GlobalEnv), add = TRUE)
  force(expr)
}

## --------------------------- the ten schemes -------------------------
## Second argument of N(0, v) in the paper is a VARIANCE, so sd = sqrt(v).
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

## ------------------------- statistic adapters ------------------------
## Tolerant accessor for packages that return a test object rather than a
## statistic.  Take the object if it is already one finite number, else the
## first candidate field holding one, else the maximum of the first
## candidate field holding a numeric vector -- the maximum because every
## statistic here is oriented so that larger means more dependence.  NA if
## nothing matches, which self_test() surfaces as a failure.
##
## To correct the candidate lists, look at what the packages actually return:
##   set.seed(1); X <- cbind(runif(100), runif(100))
##   str(BET::BET(X, d = 3)); str(BET::BEAST(X, d = 3))
##   str(HHG::hhg.univariate.ind.stat(X[,1], X[,2], variant = "ADP-EQP-ML"))
pluck1 <- function(obj, fields) {
  if (is.numeric(obj) && length(obj) == 1L && is.finite(obj)) return(as.numeric(obj))
  for (f in fields) {
    v <- tryCatch(obj[[f]], error = function(e) NULL)
    if (is.numeric(v) && length(v) == 1L && is.finite(v)) return(as.numeric(v))
  }
  for (f in fields) {
    v <- tryCatch(obj[[f]], error = function(e) NULL)
    if (is.numeric(v) && length(v) > 1L && any(is.finite(v))) return(max(v, na.rm = TRUE))
  }
  NA_real_
}

make_stats <- function() {
  s <- list(
    B_n = function(x, y, p, q) indeptest(x, y, p = p, q = q, basis = "poly",
                                         test = "Bartlett", ties = "first")$B_stat,
    P_n = function(x, y, p, q) indeptest(x, y, p = p, q = q, basis = "poly",
                                         test = "Pillai",   ties = "first")$P_stat,
    Hoeffding = function(x, y, p, q) Hmisc::hoeffd(cbind(x, y))$D[1, 2],
    dCov = function(x, y, p, q) energy::dcor(x, y),
    HHG  = function(x, y, p, q) HHG::hhg.test(as.matrix(dist(x)), as.matrix(dist(y)),
                                              nr.perm = 0L)$sum.chisq,
    ## B = 0 breaks the package (t(xic) %*% M0); the statistic is invariant
    ## to B, verified identical at B = 2, 5, 20, so the smallest workable B
    ## is used.  $stat is a list, hence $stat$Sn and not $stat["Sn"].
    Genest = function(x, y, p, q)
      MixedIndTests::TestIndCopula(cbind(x, y), B = 2)$stat$Sn
  )
  if (USE_MIC) s$MIC <- function(x, y, p, q) minerva::mine_stat(x, y, measure = "mic")

  if (USE_NEW) {
    ## BET, binary expansion testing (Zhang, JASA 2019).  The closest
    ## competitor in design: like ours it projects the copula onto a finite
    ## basis and tests the projection, but the basis is binary rather than
    ## smooth, and the aggregation is a maximum over interactions rather
    ## than a sum over a chosen span.  The comparison reads as binary versus
    ## smooth basis at matched resolution, which is the honest framing.
    s$BET <- function(x, y, p, q)
      pluck1(BET::BET(cbind(x, y), d = d_rule(length(x))),
             c("Statistic", "statistic", "Extreme.Asymmetry",
               "extreme.asymmetry", "Max.Asymmetry", "bet.s"))

    ## BEAST, the adaptive-weight version, designed for near-uniform power
    ## and the most likely of the four to beat us somewhere.  If it does,
    ## that is the result: say so, and let the calibration and cost
    ## comparison carry the argument instead.  It resamples internally to
    ## build its weights, so its value is random given the data; that extra
    ## variability belongs to the method, not to this harness, and it is why
    ## detect_rank_invariant() sets a common seed before each of its calls.
    s$BEAST <- function(x, y, p, q)
      pluck1(BET::BEAST(cbind(x, y), d = d_rule(length(x))),
             c("BEAST.Statistic", "beast.statistic", "Statistic", "statistic"))

    ## Chatterjee's xi, ASYMMETRIC by construction, so the symmetrised
    ## maximum is the fair comparison for a test of independence.  Shi,
    ## Drton and Han (Biometrika 109, 2022) show it has no power against
    ## root-n local alternatives, so it should lose at these sample sizes --
    ## worth reporting, being a consequence of design rather than tuning,
    ## and the same kind of statement Proposition 1 makes about our own
    ## blind directions.
    s$xi <- function(x, y, p, q) {
      f <- function(a, b) pluck1(XICOR::xicor(a, b), c("xi", "coef", "statistic"))
      max(f(x, y), f(y, x))
    }

    ## Fast HHG: NOT a faster algorithm for the HHG statistic above but a
    ## DIFFERENT statistic, splitting only at nr.atoms equidistant order
    ## statistics instead of at every observation, giving O(nr.atoms^4),
    ## constant in n once the data are ranked.  That is why it scales, and
    ## the price is a coarsened partition.  Report it ALONGSIDE the original,
    ## never instead: at n = 50 with 40 atoms on 50 points the two are nearly
    ## the same test, by n = 200 they are not, and the gap is informative.
    ## The package recommends n > 100, so the n = 50 column is outside its
    ## intended regime and must be flagged as such in the paper.
    s$HHG_fast <- function(x, y, p, q) {
      na <- min(40L, max(10L, length(x) %/% 2L))
      pluck1(HHG::hhg.univariate.ind.stat(x, y, variant = "ADP-EQP-ML",
                                          nr.atoms = na),
             c("statistic", "MinP", "sum.chisq", "m.stats"))
    }
  }
  s
}

## Does the statistic depend on the data only through the ranks?  If so it
## is distribution-free under independence and needs calibrating once per
## (n, p, q) rather than once per cell.  Checked, not assumed.  A common
## seed is set before each of the two calls so that a statistic with
## internal randomness, such as BEAST, is compared fairly rather than being
## rejected for its own noise.
detect_rank_invariant <- function(n = 100L, p = 3L, q = 3L, reps = 3L, tol = 1e-8) {
  S <- make_stats(); ok <- rep(TRUE, length(S)); names(ok) <- names(S)
  for (r in seq_len(reps)) {
    set.seed(1000L + r)
    x <- rnorm(n); y <- 0.6*x + rnorm(n)
    rx <- rank(x); ry <- rank(y)          # ranks are a monotone transform
    for (nm in names(S)) {
      a <- tryCatch({ set.seed(99L); as.numeric(S[[nm]](x,  y,  p, q)) },
                    error = function(e) NA_real_)
      b <- tryCatch({ set.seed(99L); as.numeric(S[[nm]](rx, ry, p, q)) },
                    error = function(e) NA_real_)
      if (!isTRUE(abs(a - b) <= tol * max(1, abs(a)))) ok[nm] <- FALSE
    }
  }
  ok
}

self_test <- function(n, p, q) {
  set.seed(1); x <- runif(n); y <- x + 0.5*rnorm(n)
  S <- make_stats(); ok <- logical(length(S)); names(ok) <- names(S)
  secs <- setNames(rep(NA_real_, length(S)), names(S))
  for (nm in names(S)) {
    t0 <- Sys.time()
    v <- tryCatch(as.numeric(S[[nm]](x, y, p, q)), error = function(e) e)
    secs[nm] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok[nm] <- !(inherits(v, "error") || length(v) != 1L || !is.finite(v))
    message(sprintf("  %-10s %s", nm, if (ok[nm]) sprintf("ok  %-12.6g (%.4f s)", v, secs[nm])
                    else paste("FAILED:", if (inherits(v, "error")) conditionMessage(v) else "bad value")))
  }
  attr(ok, "secs") <- secs
  ok
}

## ------------------------- startup diagnostics -----------------------
message("adapter self-test at n = ", N_GRID[1], ":")
STAT_OK <- self_test(N_GRID[1], p_rule(N_GRID[1]), p_rule(N_GRID[1]))
SECS    <- attr(STAT_OK, "secs")
if (!all(STAT_OK)) warning("dropping: ", paste(names(STAT_OK)[!STAT_OK], collapse = ", "))
KEEP <- names(STAT_OK)[STAT_OK]
SECS <- SECS[KEEP]
stopifnot("B_n must be available" = "B_n" %in% KEEP)

message("\nrank-invariance detection (distribution-free under H0):")
RI <- detect_rank_invariant()[KEEP]
for (nm in KEEP) message(sprintf("  %-10s %s", nm,
  if (RI[nm]) "rank-invariant: calibrated once per (n, p, q)" else "value-based: calibrated per cell"))
DF_NAMES  <- names(RI)[RI]; CELL_NAMES <- names(RI)[!RI]

## ------------------------- cost projection ---------------------------
PQ <- unique(do.call(rbind, lapply(N_GRID, function(n)
  do.call(rbind, lapply(unique(c(p_rule(n), if (SECONDARY_Q1) 1L)),
                        function(q) c(n, p_rule(n), q))))))
NCELL <- length(N_GRID) * length(gen) * length(NU) *
  (1L + as.integer(SECONDARY_Q1) * length(FUNCTIONAL) / length(gen))
proj <- data.frame(test = KEEP, secs = as.numeric(SECS[KEEP]), stringsAsFactors = FALSE)
proj$evals <- vapply(KEEP, function(t)
  NCELL * cap_for(t, R1) +
    if (RI[t]) nrow(PQ) * cap_for(t, R0_DF) else NCELL * cap_for(t, R0), 0)
proj$hours <- proj$secs * proj$evals / 3600
proj <- proj[order(-proj$hours), ]
message("\nprojected cost, from the self-test timings:")
print(proj, row.names = FALSE, digits = 3)
message(sprintf("total %.1f single-core hours, %.1f h on %d cores",
                sum(proj$hours), sum(proj$hours)/NCORES, NCORES))
if (sum(proj$hours) > MAX_PROJECTED_HOURS)
  stop(sprintf(paste("projected %.0f single-core hours exceeds MAX_PROJECTED_HOURS = %.0f.",
                     "Set EVAL_CAP for the statistics at the top of that table, or raise",
                     "the limit deliberately."), sum(proj$hours), MAX_PROJECTED_HOURS))

## ------------- critical values for the distribution-free ones --------
## A cached table is reused only if it covers every statistic now kept;
## tables written before the new competitors were added do not, and are
## recomputed rather than silently leaving gaps.
crit_df_table <- list()
for (r in seq_len(nrow(PQ))) {
  n <- PQ[r, 1]; p <- PQ[r, 2]; q <- PQ[r, 3]
  key <- sprintf("n%d_p%d_q%d", n, p, q)
  f <- file.path(OUTDIR, paste0("critdf_", key, ".rds"))
  if (file.exists(f)) {
    cv <- tryCatch(readRDS(f), error = function(e) NULL)
    if (!is.null(cv) && all(DF_NAMES %in% names(cv))) { crit_df_table[[key]] <- cv; next }
    message(sprintf("  stale critical-value table for %s (missing %s): recomputing", key,
                    paste(setdiff(DF_NAMES, names(cv)), collapse = ", ")))
  }
  message(sprintf("calibrating distribution-free statistics, %s, R = %g", key, R0_DF))
  S <- make_stats()[DF_NAMES]
  set.seed(SEED + 31L*n + 7L*p + q)
  caps <- vapply(DF_NAMES, cap_for, 0L, default = R0_DF)
  M <- matrix(NA_real_, max(caps), length(S), dimnames = list(NULL, DF_NAMES))
  for (i in seq_len(max(caps))) {
    x <- runif(n); y <- runif(n)             # independent; marginals irrelevant
    with_frozen_rng(
      for (k in DF_NAMES) if (i <= caps[k]) M[i, k] <- as.numeric(S[[k]](x, y, p, q)))
  }
  cv <- vapply(DF_NAMES, function(k)
    quantile(M[seq_len(caps[k]), k], 1 - ALPHA, names = FALSE, na.rm = TRUE), 0)
  saveRDS(cv, f); crit_df_table[[key]] <- cv
}

## ------------------------------ one cell -----------------------------
cell_file <- function(cell) file.path(OUTDIR,
  sprintf("cell_n%d_s%02d_nu%03d_q%d.rds", cell$n, cell$si, round(cell$nu*100), cell$q))

## A cached cell counts as done only if it holds every statistic now kept.
cell_done <- function(cell, keep) {
  f <- cell_file(cell)
  if (!file.exists(f)) return(FALSE)
  d <- tryCatch(readRDS(f), error = function(e) NULL)
  !is.null(d) && all(keep %in% d$test)
}

run_cell <- function(cell) {
  if (cell_done(cell, KEEP)) return(invisible(NULL))
  library(dependence)
  S <- make_stats()[KEEP]
  n <- cell$n; nu <- cell$nu; p <- cell$p; q <- cell$q; g <- gen[[cell$scheme]]
  set.seed(SEED + 1009L*n + 7919L*cell$si + 104729L*as.integer(round(nu*10)) + 13L*q)

  crit <- crit_df_table[[sprintf("n%d_p%d_q%d", n, p, q)]]
  if (length(CELL_NAMES)) {                 # per-cell calibration, same marginals
    caps <- vapply(CELL_NAMES, cap_for, 0L, default = R0)
    null <- matrix(NA_real_, max(caps), length(CELL_NAMES),
                   dimnames = list(NULL, CELL_NAMES))
    for (i in seq_len(max(caps))) {
      d <- g(n, nu); d$y <- sample(d$y)
      with_frozen_rng(
        for (k in CELL_NAMES) if (i <= caps[k])
          null[i, k] <- as.numeric(S[[k]](d$x, d$y, p, q)))
    }
    crit <- c(crit, vapply(CELL_NAMES, function(k)
      quantile(null[seq_len(caps[k]), k], 1 - ALPHA, names = FALSE, na.rm = TRUE), 0))
  }
  crit <- crit[KEEP]

  caps1 <- vapply(KEEP, cap_for, 0L, default = R1)
  alt <- matrix(NA_real_, max(caps1), length(KEEP), dimnames = list(NULL, KEEP))
  for (i in seq_len(max(caps1))) {
    d <- g(n, nu)
    with_frozen_rng(
      for (k in KEEP) if (i <= caps1[k]) alt[i, k] <- as.numeric(S[[k]](d$x, d$y, p, q)))
  }
  pw <- vapply(KEEP, function(k)
    mean(alt[seq_len(caps1[k]), k] > crit[k], na.rm = TRUE), 0)

  out <- data.frame(n = n, scheme = cell$scheme, nu = nu, p = p, q = q,
                    test = KEEP, power = as.numeric(pw),
                    se = sqrt(as.numeric(pw)*(1 - as.numeric(pw))/caps1[KEEP]),
                    crit = as.numeric(crit),
                    R0 = ifelse(KEEP %in% DF_NAMES, R0_DF, R0), R1 = caps1[KEEP],
                    df_calibrated = KEEP %in% DF_NAMES,
                    row.names = NULL, stringsAsFactors = FALSE)
  saveRDS(out, cell_file(cell)); invisible(NULL)
}

## ------------------------------- driver ------------------------------
cells <- list()
for (n in N_GRID) for (si in seq_along(gen)) for (nu in NU) {
  sch <- names(gen)[si]
  cells[[length(cells)+1L]] <- list(n=n, si=si, scheme=sch, nu=nu, p=p_rule(n), q=p_rule(n))
  if (SECONDARY_Q1 && sch %in% FUNCTIONAL)
    cells[[length(cells)+1L]] <- list(n=n, si=si, scheme=sch, nu=nu, p=p_rule(n), q=1L)
}
cells <- cells[order(-vapply(cells, function(z) z$n, 0))]   # expensive first
done  <- vapply(cells, cell_done, TRUE, keep = KEEP)
message(sprintf("\n%d cells, %d already done with the current statistics, %d to run, %d cores",
                length(cells), sum(done), sum(!done), NCORES))
if (sum(done) == 0 && length(list.files(OUTDIR, "^cell_.*rds$")))
  message("  (existing cells predate the new competitors and will be recomputed)")

todo <- cells[!done]
if (length(todo)) {
  t0 <- Sys.time(); run_cell(todo[[1]])
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  message(sprintf("first cell (n = %d) took %.0f s; remaining %d cells approx %.1f h on %d cores",
                  todo[[1]]$n, dt, length(todo)-1, dt*(length(todo)-1)/3600/NCORES, NCORES))
  cl <- makeCluster(NCORES)
  clusterExport(cl, c("gen","make_stats","run_cell","cell_file","cell_done","pluck1",
                      "with_frozen_rng","N_GRID","R0","R0_DF","R1","ALPHA","SEED",
                      "USE_MIC","USE_NEW","p_rule","d_rule","cap_for","EVAL_CAP",
                      "KEEP","DF_NAMES","CELL_NAMES","crit_df_table","OUTDIR"))
  invisible(clusterEvalQ(cl, library(dependence)))
  invisible(parLapplyLB(cl, todo[-1], run_cell))
  stopCluster(cl)
}

## ----------------------------- assemble ------------------------------
res <- do.call(rbind, lapply(list.files(OUTDIR, "^cell_.*rds$", full.names = TRUE), readRDS))
write.csv(res, paste0(OUTSTEM, ".csv"), row.names = FALSE)
prim <- res[res$q == res$p, ]

## The distribution-free statistics now take one critical value per
## (n, p, q) by construction, so this check is only meaningful for the ones
## calibrated per cell, where a large spread signals an unstable statistic.
for (tst in intersect(CELL_NAMES, unique(prim$test))) for (n in N_GRID) {
  v <- prim$crit[prim$test == tst & prim$n == n]
  if (length(v) && all(is.finite(v)))
    message(sprintf("  per-cell critical values: %-10s n=%-4d spread %.1e", tst, n,
                    diff(range(v))/mean(v)))
}
summ <- do.call(rbind, lapply(split(prim, list(prim$n, prim$test), drop = TRUE),
  function(d) data.frame(n = d$n[1], test = d$test[1],
                         mean_power = mean(d$power), sd_power = sd(d$power))))
rk <- do.call(rbind, lapply(split(prim, list(prim$n, prim$scheme, prim$nu), drop = TRUE),
  function(d) data.frame(n = d$n[1], test = d$test, rank = rank(-d$power, ties.method = "average"))))
summ$mean_rank <- vapply(seq_len(nrow(summ)), function(i)
  mean(rk$rank[rk$n == summ$n[i] & rk$test == summ$test[i]]), 0)
summ <- summ[order(summ$n, summ$mean_rank), ]
print(summ, row.names = FALSE)
write.csv(summ, paste0(OUTSTEM, "_summary.csv"), row.names = FALSE)

## LaTeX fragment for the main paper
fmt <- function(v) sub("0\\.", "0{\\\\cdot}", sprintf("%.3f", v))
TESTS <- unique(summ$test[summ$n == N_GRID[1]])
con <- file(paste0(OUTSTEM, "_summary.tex"), open = "wt")
writeLines(c("\\begin{table}",
"\\tbl{Size-adjusted empirical power across ten dependence schemes and twenty noise levels}{%",
"\\begin{tabular}{l rrr rrr rrr}", "\\hline",
paste0(" & ", paste(sprintf("\\multicolumn{3}{c}{$n = %d$}", N_GRID), collapse = " & "), " \\\\"),
paste0("Test & ", paste(rep("Mean & SD & Rank", length(N_GRID)), collapse = " & "), " \\\\"),
"\\hline",
vapply(TESTS, function(t) {
  cc <- vapply(N_GRID, function(n) { d <- summ[summ$n == n & summ$test == t, ]
    if (!nrow(d)) return("--- & --- & ---")
    sprintf("$%s$ & $%s$ & $%s$", fmt(d$mean_power), fmt(d$sd_power), fmt(d$mean_rank)) }, "")
  sprintf("%s & %s \\\\", gsub("_", "\\\\_", t), paste(cc, collapse = " & ")) }, ""),
"\\hline", "\\end{tabular}}", "\\begin{tabnote}",
sprintf(paste("Critical values are size-adjusted from $%g$ null datasets per cell, obtained by permuting $Y$,",
              "except for the statistics depending on the data only through the ranks, which are",
              "distribution-free under independence and are calibrated once per $(n,p,q)$ from $%g$ null",
              "datasets. Power is over $%g$ independent draws from the alternative, so the largest standard",
              "error is $%s$. Orders follow rule (12), giving $p=q=%s$; the binary expansion depth for BET and",
              "BEAST is $d=%s$, matched to the number of directions examined. The mean rank is the average",
              "rank, one being best, within each combination."),
        R0, R0_DF, R1, fmt(0.5/sqrt(R1)),
        paste(vapply(N_GRID, p_rule, 0L), collapse = ", "),
        paste(vapply(N_GRID, d_rule, 0L), collapse = ", ")),
"\\end{tabnote}", "\\label{tab:power_summary}", "\\end{table}"), con)
close(con)
message("\nWritten: ", OUTSTEM, ".csv, ", OUTSTEM, "_summary.csv, ", OUTSTEM, "_summary.tex")

## =====================================================================
##  The second table, which is the one that answers the referee
##
##  Power is not the axis on which this test wins, and a power table alone
##  is what drew the objection that the advance over existing work was not
##  established.  The ordering is not close on calibration and cost, so
##  record those: seconds per test as n grows, and whether a calibrated
##  p-value is available with no resampling at run time.  Of the statistics
##  here only B_n and P_n have a closed-form reference distribution; BET has
##  an exact null, the fast HHG has precomputed null tables that must
##  themselves be simulated per n, and BEAST, HHG, dCov, Genest and MIC have
##  nothing.  Each statistic stops being timed once it exceeds a minute.
## =====================================================================
RUN_SCALING <- TRUE
if (RUN_SCALING) {
  S <- make_stats()[KEEP]
  ns <- c(1e3, 1e4, 1e5, 1e6)
  sc <- expand.grid(test = KEEP, n = ns, secs = NA_real_, stringsAsFactors = FALSE)
  for (k in seq_len(nrow(sc))) {
    n <- as.integer(sc$n[k]); t <- sc$test[k]
    prev <- sc$secs[sc$test == t & sc$n < n]
    if (length(prev) && any(is.na(prev) | prev > 60)) next
    set.seed(11); x <- rnorm(n); y <- 0.5*x + rnorm(n); p <- p_rule(n)
    t0 <- Sys.time()
    okk <- tryCatch({ S[[t]](x, y, p, p); TRUE }, error = function(e) FALSE)
    sc$secs[k] <- if (okk) as.numeric(difftime(Sys.time(), t0, units = "secs")) else NA_real_
  }
  print(sc[order(sc$test, sc$n), ], row.names = FALSE)
  write.csv(sc, paste0(OUTSTEM, "_scaling.csv"), row.names = FALSE)
  message("Written: ", OUTSTEM, "_scaling.csv")
}
