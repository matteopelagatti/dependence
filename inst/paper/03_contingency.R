## =====================================================================
##  03_contingency.R
##
##  Size-adjusted power on two-way tables whose categories are genuinely
##  NOMINAL.  Replaces the study behind Figure 3 of the paper.
##
##  WHY.  The original design sampled a bivariate normal and discretized the
##  margins, so the categories arrived in their natural order and every rank-
##  or distance-based competitor could exploit that ordering, while the
##  indicator basis -- which is Pearson's chi-squared, by Theorem 6 -- could
##  not.  That is why Hoeffding overtook our test on 5 x 5 tables at
##  rho >= 0.2.  Here the labels are permuted by a fixed arbitrary
##  permutation after discretization, so the ordering carries no information
##  and every method faces the same problem.
##  Permutations: 3 x 3 -> (2, 3, 1);  5 x 5 -> (3, 5, 2, 1, 4).
##
##  COMPETITORS.  Two of them, Hoeffding and the Genest copula functional,
##  work from ranks; two, distance covariance and HHG, work from distances on
##  the integer labels.  All four are therefore sensitive to the relabelling,
##  which is the point: the figure is meant to show what their advantage on
##  the original ordinal design was worth.  Pearson's chi-squared is the
##  invariant reference, and coincides with P_n wherever rule (12) saturates.
##
##  SATURATION.  Rule (12) reaches the c-1 ceiling of section 3.1 in half the
##  cells: 3 x 3 gives p = 1, 2, 2, 2 and 5 x 5 gives p = 1, 2, 3, 4 at
##  n = 25, 50, 100, 200.  Where it saturates, the polynomial expansion of
##  the numeric labels spans exactly the indicator space, so by Remark 5 the
##  two constructions give the identical statistic and the relabelling can
##  have no effect at all.  There is then nothing to be gained by computing
##  both, so B_n and P_n are obtained from the indicator basis in those cells
##  and from the polynomial basis on the labels only where p < c-1.
##
##  TIES.  Remark 5 requires tied observations to SHARE a rank; any
##  convention that separates them lets the within-category rank spread enter
##  as noise and breaks the equivalence.  This script therefore uses
##  ties = "average".  With mid-ranks the scaled ranks take only c distinct
##  values, so p must not exceed c-1 or the design is exactly singular; the
##  cap is enforced, not optional.
##
##  The equivalence is verified numerically at startup rather than assumed:
##  see check_remark5(), which also re-derives Theorem 6.
##
##  Size adjustment is cheap: the marginals of a discretized normal do not
##  depend on rho, so one null calibration per (table size, n) serves every
##  rho.  The rho = 0 column is then a check on the calibration, not a
##  result -- it must come out at the nominal level by construction.
## =====================================================================

library(dependence)
library(parallel)

## ---------------------------- configuration --------------------------
DIMS    <- list("3x3" = list(k = 3L, cut = c(-0.6, 0.6),            perm = c(2L, 3L, 1L)),
                "5x5" = list(k = 5L, cut = c(-1.0, -0.3, 0.3, 1.0), perm = c(3L, 5L, 2L, 1L, 4L)))
N_GRID  <- c(25L, 50L, 100L, 200L)
RHO     <- seq(0, 0.9, by = 0.1)
ALPHA   <- 0.05
R0      <- 2e4          # null datasets per (table size, n)
R1      <- 1e4          # datasets per (table size, n, rho)
INCLUDE_CHISQ <- TRUE   # relabelling-invariant reference; dropped from the
# figure automatically wherever it duplicates P_n
TIES    <- "average"    # mandatory for Remark 5; see the header
SEED    <- 20260811L
NCORES  <- max(1L, detectCores() - 1L)
OUTDIR  <- "inst/paper/results"
OUTSTEM <- "contingency_nominal"
dir.create(OUTDIR, showWarnings = FALSE)

p_rule <- function(n) max(1L, as.integer(floor(n^0.31)) - 1L)
## With mid-ranks the scaled ranks take only k distinct values, so p > k-1
## makes S_uu exactly singular: the cap is a requirement, not a preference.
orders <- function(n, k) min(p_rule(n), k - 1L)
saturated <- function(n, k) orders(n, k) >= k - 1L

## ------------------------------ the data -----------------------------
gen_table <- function(n, rho, spec) {
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  list(x = spec$perm[findInterval(z1, spec$cut) + 1L],
       y = spec$perm[findInterval(z2, spec$cut) + 1L])
}

## --------------------- Remark 5 and Theorem 6, checked ----------------
check_remark5 <- function(n = 600L, tol = 1e-8) {
  for (dn in names(DIMS)) {
    spec <- DIMS[[dn]]; k <- spec$k
    set.seed(7); d <- gen_table(n, 0.4, spec)
    a  <- indeptest(d$x, d$y, p = k - 1L, q = k - 1L, basis = "poly",
                    test = c("Pillai", "Bartlett"), ties = "average")
    b  <- indeptest(factor(d$x), factor(d$y), basis = "dummy",
                    test = c("Pillai", "Bartlett"))
    cs <- suppressWarnings(as.numeric(
      chisq.test(table(d$x, d$y), correct = FALSE)$statistic))
    dP <- abs(a$P_stat - b$P_stat); dB <- abs(a$B_stat - b$B_stat)
    dC <- abs(b$P_stat - cs)
    message(sprintf("  %-4s |P_poly - P_dummy| = %.2e   |B_poly - B_dummy| = %.2e   |P_dummy - X2| = %.2e",
                    dn, dP, dB, dC))
    if (max(dP, dB, dC) > tol * max(1, abs(cs)))
      stop("Remark 5 / Theorem 6 does not hold numerically at p = c-1 with ",
           "ties = 'average'.  Either the tie handling or the basis is wrong; ",
           "fix that before running the study.")
  }
  invisible(TRUE)
}

## ------------------------- statistic adapters ------------------------
## Each takes the integer labels, the order p and the number of levels k,
## and returns one number, larger meaning more dependence.
make_stats <- function() {
  s <- list(
    B_n = function(x, y, p, k)
      if (p >= k - 1L)                       # saturated: identical anyway
        indeptest(factor(x), factor(y), basis = "dummy", test = "Bartlett")$B_stat
    else
      indeptest(x, y, p = p, q = p, basis = "poly", test = "Bartlett",
                ties = TIES)$B_stat,
    P_n = function(x, y, p, k)
      if (p >= k - 1L)
        indeptest(factor(x), factor(y), basis = "dummy", test = "Pillai")$P_stat
    else
      indeptest(x, y, p = p, q = p, basis = "poly", test = "Pillai",
                ties = TIES)$P_stat,
    Hoeffding = function(x, y, p, k) Hmisc::hoeffd(cbind(x, y))$D[1, 2],
    ## dCov and HHG are here precisely BECAUSE they are not invariant to a
    ## relabelling of the categories: both read the integer codes as points on
    ## the line, so that under the 5 x 5 permutation categories 3 and 5 are
    ## held four apart and 3 and 2 one apart, when neither statement means
    ## anything.  That sensitivity is the object of this study, so excluding
    ## them would remove the comparison the redesign was built to make.
    dCov = function(x, y, p, k) energy::dcor(x, y),
    HHG  = function(x, y, p, k) HHG::hhg.test(as.matrix(dist(x)), as.matrix(dist(y)),
                                              nr.perm = 0L)$sum.chisq,
    ## Two statistics from the continuous study are omitted, deliberately.
    ## MIC searches over grid partitions of the plane and is degenerate with
    ## three to five distinct values.  The Bakirov, Rizzo and Szekely (2006)
    ## coefficient, exposed by energy::indep.test(method = "mvI"), is a
    ## fourfold sum: timed at 3504 ms per dataset at n = 200 against some 50 ms
    ## for everything retained here together, it alone accounted for 99.7 per
    ## cent of the run, forty hours against two, and dCov already represents
    ## the energy family.
    Genest = function(x, y, p, k)
      MixedIndTests::TestIndCopula(cbind(x, y), B = 2)$stat$Sn
  )
  if (INCLUDE_CHISQ)
    s$ChiSq <- function(x, y, p, k)
      suppressWarnings(as.numeric(chisq.test(table(x, y), correct = FALSE)$statistic))
  s
}

self_test <- function() {
  set.seed(1); spec <- DIMS[["5x5"]]; d <- gen_table(200L, 0.5, spec)
  p <- orders(200L, spec$k); S <- make_stats()
  ok <- logical(length(S)); names(ok) <- names(S)
  for (nm in names(S)) {
    t0 <- Sys.time()
    v <- tryCatch(as.numeric(S[[nm]](d$x, d$y, p, spec$k)), error = function(e) e)
    dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ok[nm] <- !(inherits(v, "error") || length(v) != 1L || !is.finite(v))
    message(sprintf("  %-10s %s", nm, if (ok[nm]) sprintf("ok  %-12.6g (%.4f s)", v, dt)
                    else paste("FAILED:", if (inherits(v, "error")) conditionMessage(v) else "bad value")))
  }
  for (nm in intersect("Genest", names(S)[ok])) {
    v <- vapply(1:3, function(i) as.numeric(S[[nm]](d$x, d$y, p, spec$k)), 0)
    if (diff(range(v)) > 1e-8)
      warning(nm, " depends on its internal resampling; the extraction is wrong")
  }
  ok
}

## ------------------------- calibration and cells ---------------------
calibrate <- function(dim_nm, n) {
  f <- file.path(OUTDIR, sprintf("crit_%s_n%d.rds", dim_nm, n))
  if (file.exists(f)) return(readRDS(f))
  spec <- DIMS[[dim_nm]]; p <- orders(n, spec$k); S <- make_stats()[KEEP]
  set.seed(SEED + 101L * n + 7L * spec$k)
  M <- matrix(NA_real_, R0, length(S), dimnames = list(NULL, names(S)))
  for (i in seq_len(R0)) { d <- gen_table(n, 0, spec)
  M[i, ] <- vapply(S, function(g) as.numeric(g(d$x, d$y, p, spec$k)), 0) }
  cv <- apply(M, 2L, quantile, probs = 1 - ALPHA, names = FALSE, na.rm = TRUE)
  saveRDS(cv, f); cv
}

run_cell <- function(cell) {
  f <- file.path(OUTDIR, sprintf("cell_%s_n%d_rho%02d.rds", cell$dim, cell$n,
                                 round(cell$rho * 10)))
  if (file.exists(f)) return(invisible(NULL))
  library(dependence)
  spec <- DIMS[[cell$dim]]; k <- spec$k; p <- orders(cell$n, k)
  S <- make_stats()[KEEP]; crit <- CRIT[[paste(cell$dim, cell$n)]]
  set.seed(SEED + 1009L*cell$n + 7919L*k + 104729L*as.integer(round(cell$rho*10)))
  A <- matrix(NA_real_, R1, length(S), dimnames = list(NULL, names(S)))
  for (i in seq_len(R1)) { d <- gen_table(cell$n, cell$rho, spec)
  A[i, ] <- vapply(S, function(g) as.numeric(g(d$x, d$y, p, k)), 0) }
  pw <- colMeans(sweep(A, 2L, crit[names(S)], ">"), na.rm = TRUE)
  out <- data.frame(dim = cell$dim, n = cell$n, rho = cell$rho, p = p,
                    saturated = saturated(cell$n, k), test = names(S),
                    power = as.numeric(pw),
                    se = sqrt(as.numeric(pw)*(1 - as.numeric(pw))/R1),
                    crit = as.numeric(crit[names(S)]), R0 = R0, R1 = R1,
                    row.names = NULL, stringsAsFactors = FALSE)
  saveRDS(out, f); invisible(NULL)
}

## ------------------------------- driver ------------------------------
message("Remark 5 and Theorem 6, verified numerically at p = q = c-1:")
check_remark5()

message("\nadapter self-test:")
STAT_OK <- self_test()
if (!all(STAT_OK)) warning("dropping: ", paste(names(STAT_OK)[!STAT_OK], collapse = ", "))
KEEP <- names(STAT_OK)[STAT_OK]
stopifnot("B_n must be available" = "B_n" %in% KEEP)

message("\norders in use, and which construction B_n and P_n take:")
for (dn in names(DIMS)) { k <- DIMS[[dn]]$k
message(sprintf("  %-4s %s", dn, paste(sprintf("n=%d: p=%d %s", N_GRID,
                                               vapply(N_GRID, orders, 0L, k = k),
                                               ifelse(vapply(N_GRID, saturated, TRUE, k = k), "[indicator]", "[polynomial]")),
                                       collapse = "   "))) }

CRIT <- list()
for (dn in names(DIMS)) for (n in N_GRID) {
  message(sprintf("calibrating %s, n = %d, R0 = %g", dn, n, R0))
  CRIT[[paste(dn, n)]] <- calibrate(dn, n)
}

cells <- list()
for (dn in names(DIMS)) for (n in N_GRID) for (rho in RHO)
  cells[[length(cells) + 1L]] <- list(dim = dn, n = n, rho = rho)
cells <- cells[order(-vapply(cells, function(z) z$n, 0))]
message(sprintf("\n%d cells, %d cores", length(cells), NCORES))

cl <- makeCluster(NCORES)
clusterExport(cl, c("DIMS","N_GRID","RHO","ALPHA","R0","R1","SEED","TIES","KEEP","CRIT",
                    "OUTDIR","INCLUDE_CHISQ","p_rule","orders","saturated",
                    "gen_table","make_stats","run_cell"))
invisible(clusterEvalQ(cl, library(dependence)))
invisible(parLapplyLB(cl, cells, run_cell))
stopCluster(cl)

## ----------------------------- assemble ------------------------------
res <- do.call(rbind, lapply(list.files(OUTDIR, "^cell_.*rds$", full.names = TRUE), readRDS))
res <- res[order(res$dim, res$n, res$rho, res$test), ]
write.csv(res, paste0(OUTSTEM, ".csv"), row.names = FALSE)

sz <- res[res$rho == 0, ]
message("\nempirical size at rho = 0 (target ", ALPHA, ", holds by construction):")
for (dn in names(DIMS)) for (n in N_GRID) {
  s <- sz[sz$dim == dn & sz$n == n, ]
  message(sprintf("  %-4s n=%-4d %s", dn, n,
                  paste(sprintf("%s %.3f", s$test, s$power), collapse = "  ")))
}

## Which tests coincide, cell by cell?  Powers are exact means, so a
## duplicate pair differs by exactly zero.
wide <- reshape(res[, c("dim","n","rho","test","power")],
                idvar = c("dim","n","rho"), timevar = "test", direction = "wide")
tests <- sub("^power\\.", "", grep("^power\\.", names(wide), value = TRUE))
message("\nmaximum |difference| in power between pairs of tests:")
dup <- character(0)
for (i in seq_along(tests)) for (j in seq_len(i - 1L)) {
  dmax <- max(abs(wide[[paste0("power.", tests[i])]] - wide[[paste0("power.", tests[j])]]))
  if (dmax < 0.01)
    message(sprintf("  %-10s vs %-10s  %.4f%s", tests[i], tests[j], dmax,
                    if (dmax == 0) "   IDENTICAL" else ""))
  if (dmax == 0) dup <- union(dup, tests[i])
}
if (length(dup)) message("  -> dropping from the figure: ", paste(dup, collapse = ", "))

## ------------------------------- figure ------------------------------
PLOT <- setdiff(tests, dup)
LEGEND <- PLOT
LEGEND[LEGEND == "B_n"] <- "B(p)"
LEGEND[LEGEND == "P_n"] <- "P(p)"
LEGEND[LEGEND == "ChiSq"] <- "ChiSq = P(c-1)"

pdf(paste0(OUTDIR, "/", OUTSTEM, ".pdf"), width = 9.5, height = 5)
op <- par(mfrow = c(length(DIMS), length(N_GRID)), mar = c(3.2, 3.2, 2, 0.6),
          mgp = c(2, 0.6, 0), cex = 0.75)
for (dn in names(DIMS)) for (n in N_GRID) {
  plot(NA, xlim = range(RHO), ylim = c(0, 1), xlab = expression(rho),
       ylab = "rejection rate",
       main = sprintf("%s, n = %d, p = %d%s", dn, n, orders(n, DIMS[[dn]]$k),
                      if (saturated(n, DIMS[[dn]]$k)) "" else "*"))
  abline(h = ALPHA, lty = 3, col = "grey50")
  for (j in seq_along(PLOT)) {
    s <- res[res$dim == dn & res$n == n & res$test == PLOT[j], ]
    lines(s$rho, s$power, col = j, lty = j, lwd = 1.6)
  }
  if (dn == names(DIMS)[1] && n == N_GRID[1])
    legend("topleft", LEGEND, col = seq_along(PLOT), lty = seq_along(PLOT),
           lwd = 1.6, bty = "n", cex = 0.85)
}
par(op); dev.off()
message("\nAn asterisk in a panel title marks p < c-1, the only cells in which ",
        "the labelling can affect B_n and P_n.")
message("Written: ", OUTSTEM, ".csv and ", OUTSTEM, ".pdf")
