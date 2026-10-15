## =====================================================================
##  05_scan_K_calibration.R
##
##  CALIBRATION OF K(n), THE CEILING OF THE ORDER SCAN.
##
##  Section 2.6 of the paper fixes K two orders above rule (12), giving
##  K = 4, 5, 6 at n = 50, 100, 200, and says plainly that this is a default
##  fixed in advance rather than a calibrated rate.  This script is the
##  study that replaces that sentence, or confirms it.
##
##  WHY RULE (12)'s RECIPE DOES NOT TRANSFER.  Rule (12) was found by
##  locating a SIZE boundary: the largest order at which the chi-squared
##  approximation still holds to within a tenth of the nominal level.  The
##  scan has no such boundary.  Its level is exact at every K, because
##  c_alpha is read from the simulated joint null, so size carries no
##  information about K at all.  The criterion has to be a power criterion,
##  and a power criterion depends on the alternative in a way the size
##  criterion did not.  That is the whole difficulty, and it shapes
##  everything below.
##
##  WHAT THE STUDY DOES.  It separates the question into a null side, which
##  is alternative-free and brackets K from above, and a power side, which
##  chooses within the bracket.
##
##  NULL SIDE, two measurements, neither of which needs an alternative:
##
##   (a) THE DILUTION.  Under the null the profile behaves as a sum of
##       k^2 independent chi-squared(1) variables, so
##             corr{P_n(j,j), P_n(k,k)} = j/k,
##       the profile is an Ornstein-Uhlenbeck process in log time, and the
##       effective multiplicity alpha/c_alpha grows like log K rather than
##       like K.  Exploratory runs give alpha/c_alpha = 2.1 log K to within
##       0.15 over n in {100, 200, 800} and K from 3 to 28.  If that holds
##       across the grid, raising K is nearly free and the ceiling is not
##       where the level is spent.
##
##   (b) THE INFORMATION BOUNDARY.  The mean of P_n(k,k) tracks k^2 to
##       three digits at every order, but its variance leaves 2k^2 well
##       before k^2/n reaches 1: at n = 400 the ratio var/(2k^2) is 0.97 at
##       k = 6, 1.18 at k = 16 and 1.38 at k = 20.  High orders are noisier
##       than their nominal degrees of freedom, so their contribution to
##       power decays faster than the chi-squared calculation implies.
##       This is the real brake on K, it is the analogue of the boundary
##       rule (12) sits inside, and it is measurable from null simulation
##       alone.  NOTE that it does not collapse onto k^2/n -- at k^2/n =
##       0.64 the inflation is 1.00 at n = 100 and 1.18 at n = 400 -- so
##       the scaling has to be located empirically, exactly as (12)'s was.
##       That is what the moment table written by this script is for.
##
##  POWER SIDE.  Alternatives along a FREQUENCY LADDER, sin(2 pi j x) for
##  j = 1, 2, 4, 8, because frequency is the only dimension K controls,
##  with the linear scheme at one end, where the best order is k = 1, and
##  the lattice at the other, where no rank correlation exists at all.  The
##  benchmark is the ORACLE single order: the best k for that alternative,
##  that noise level and that n, chosen with knowledge of the truth.  It is
##  not a procedure anyone can run, which is the point -- it bounds what
##  any choice of order could achieve, so
##        regret(family) = power(oracle) - power(family)
##  is the quantity to minimise, and the claim it supports is "the scan
##  costs at most x against knowing the right order".
##
##  WHAT WOULD BE THE BEST OUTCOME.  Given a log K dilution and a variance
##  boundary that bites only near k^2 of order n, the regret should be flat
##  in K over a wide range.  If it is, do NOT report a tuned rate: report
##  that the procedure is insensitive to K over that range and set K at the
##  admissibility ceiling.  "Not sensitive to this parameter" is a stronger
##  claim than a calibrated one, and unlike a power-tuned rate it cannot be
##  attacked as fitted to the design that measures it.  The script prints
##  the flat region explicitly so that this can be read off.
##
##  THREE INDEPENDENT NULL BATCHES, and the reason for each:
##    batch 1  the per-order laws, which turn P_n(k,k) into p_k and supply
##             the single-order critical values;
##    batch 2  the law of min_k p_k / w_k, which fixes c_alpha;
##    batch 3  a size check, on data used for neither.
##  Using one batch for 1 and 2 would evaluate the scan on the very draws
##  that defined its own reference distribution, and it would over-reject.
##  Batch 3 exists so that this is checked rather than asserted: every
##  family and every single order should come out at ALPHA on it.
##
##  PAIRED COMPARISONS.  Every family and every single order is evaluated
##  on the SAME alternative datasets, from one profile computed per
##  dataset.  Differences between families are therefore paired, and far
##  more precise than the marginal standard error 0.007 of a single
##  rejection rate at R1 = 5000 suggests.  This is also why R1 can be
##  modest: the study is about differences, not levels.
##
##  COST.  About 40 minutes single-core for the whole grid, dominated by
##  n = 1600, so roughly six minutes on seven cores.  Everything is cached
##  per cell and the script is resumable: interrupt it and run it again.
##
##  WRITES, all into inst/paper/results:
##    scan_K_null.csv      c_alpha, effective multiplicity, attained size
##    scan_K_moments.csv   mean/k^2, var/(2k^2), deviation of corr from j/k
##    scan_K_power.csv     power of every family and every single order
##    scan_K_summary.txt   the minimax regret table and the flat region
## =====================================================================

library(dependence)
library(parallel)

if (!exists("scan_profile", envir = asNamespace("dependence"), inherits = FALSE))
  stop("this script needs dependence:::scan_profile; install the current package",
       call. = FALSE)
sp <- dependence:::scan_profile
sb <- dependence:::scan_basis

## ------------------------------- settings ----------------------------
SEED   <- 20261015L
N_GRID <- c(25L, 50L, 100L, 200L, 400L, 800L, 1600L)
R0     <- 20000L      # datasets per null batch, three batches per n
R1     <-  5000L      # datasets per alternative cell
ALPHA  <- 0.05
NU     <- c(seq(0.2, 2.0, by = 0.2), 2.5, 3.0, 4.0)
NCORES <- max(1L, detectCores() - 1L)

## The ceiling of the GRID, not of the recommended rule: the point is to
## map past the useful range so that the boundary can be seen.  sqrt(n) is
## the scale at which k^2/n reaches 1, which (b) above says is already well
## inside the degraded region.
K_max <- function(n) max(4L, min(40L, as.integer(floor(sqrt(n)))))

p_rule <- function(n) max(1L, as.integer(floor(n^0.31)) - 1L)

## ---------------------------- where output goes ----------------------
## Identical to the block in 02_power_continuous.R.  Each script carries
## its own copy on purpose: they are run independently, and a shared
## helper file would be one more thing to keep in step.
paper_dir <- function() {
  cand <- character(0)
  a <- commandArgs(trailingOnly = FALSE)
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) cand <- c(cand, dirname(sub("^--file=", "", m[1L])))
  for (i in rev(seq_len(sys.nframe()))) {
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (is.character(of) && length(of) == 1L) { cand <- c(cand, dirname(of)); break }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (is.character(p) && length(p) == 1L && nzchar(p)) cand <- c(cand, dirname(p))
  }
  cand <- c(cand, getwd())
  is_paper <- function(d)
    basename(d) == "paper" && length(list.files(d, "^0[0-9].*\\.R$")) > 0L
  for (st in cand) {
    d <- tryCatch(normalizePath(st, winslash = "/", mustWork = TRUE),
                  error = function(e) NA_character_)
    if (is.na(d)) next
    for (up in 0:5) {
      if (is_paper(d)) return(d)
      if (file.exists(file.path(d, "DESCRIPTION")))
        for (q in c(file.path(d, "inst", "paper"), file.path(d, "paper")))
          if (dir.exists(q)) return(q)
      p <- dirname(d)
      if (identical(p, d)) break
      d <- p
    }
  }
  stop("cannot locate inst/paper from \"", getwd(), "\".", call. = FALSE)
}
PAPER   <- paper_dir()
OUTDIR  <- file.path(PAPER, "results")
OUTSTEM <- file.path(OUTDIR, "scan_K")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
message("results directory: ", OUTDIR)

## saveRDS() writes in place, so an interrupted write leaves a file that
## exists and cannot be read; every later run then finds it and dies.
save_rds <- function(object, file) {
  tmp <- paste0(file, ".part", Sys.getpid())
  saveRDS(object, tmp)
  if (!file.rename(tmp, file)) {
    ok <- file.copy(tmp, file, overwrite = TRUE); unlink(tmp)
    if (!ok) stop("could not write ", file, call. = FALSE)
  }
  invisible(file)
}
read_rds <- function(file)
  if (file.exists(file)) tryCatch(readRDS(file), error = function(e) NULL) else NULL

## The results directory is shared with two other studies whose cells are
## also called cell_*.rds.  One pattern per study, and this one matches
## nothing of theirs.
NULL_PAT <- "^scanK_null_n[0-9]+\\.rds$"
CELL_PAT <- "^scanK_cell_n[0-9]+_s[0-9]{2}_nu[0-9]{3}\\.rds$"
null_file <- function(n) file.path(OUTDIR, sprintf("scanK_null_n%d.rds", n))
cell_file <- function(n, si, nu)
  file.path(OUTDIR, sprintf("scanK_cell_n%d_s%02d_nu%03d.rds", n, si, round(nu*100)))

## Datasets must depend only on the cell, never on what else is computed in
## the same run, or a resumed run is not the run it resumes.
cell_seed <- function(n, si, nu)
  SEED + 1009L*as.integer(n) + 7919L*as.integer(si) + 104729L*as.integer(round(nu*10))
null_seed <- function(n, batch) SEED + 31L*as.integer(n) + 1000003L*as.integer(batch)

## ------------------------------ the schemes --------------------------
## Second argument of N(0, v) is a VARIANCE, matching the paper, so the
## sine schemes carry sd = 2 nu as in 02_power_continuous.R; j = 2 and
## j = 4 reproduce "sine" and "sine 8pi" of the main study exactly, which
## is a free cross-check against power_continuous.csv.
gen <- list(
  "linear"     = function(n, nu) { x <- runif(n); list(x = x, y = x + nu*rnorm(n)) },
  "sine 2pi"   = function(n, nu) { x <- runif(n); list(x = x, y = sin(2*pi*x)  + nu*rnorm(n, 0, 2)) },
  "sine 4pi"   = function(n, nu) { x <- runif(n); list(x = x, y = sin(4*pi*x)  + nu*rnorm(n, 0, 2)) },
  "sine 8pi"   = function(n, nu) { x <- runif(n); list(x = x, y = sin(8*pi*x)  + nu*rnorm(n, 0, 2)) },
  "sine 16pi"  = function(n, nu) { x <- runif(n); list(x = x, y = sin(16*pi*x) + nu*rnorm(n, 0, 2)) },
  "lattice"    = function(n, nu) {
    i <- sample.int(4L, n, replace = TRUE)
    j <- ifelse(i == 1L | i == 4L, sample(c(2L, 3L), n, TRUE), sample(c(1L, 4L), n, TRUE))
    list(x = (i - runif(n))/4, y = (j - runif(n))/4 + nu*rnorm(n, 0, 0.15))
  })

## ------------------------------ the families -------------------------
## A family is a set of orders and a set of weights summing to one; it
## rejects when min_k p_k / w_k <= c_alpha.  Uniform weights over a dense
## block is the version in the paper; the geometric ladders exploit the
## fact that corr(j, k) = j/k, so orders 1, 2, 4, 8, ... are equally spaced
## in log time and reach the same top order far more cheaply; the weighted
## families use every order up to the ceiling with summable weights, which
## removes the need to choose K at all and is the design most worth
## settling.  The single order of rule (12) is carried as the reference.
families <- function(Km, n) {
  f <- list()
  add <- function(name, k, w = NULL) {
    k <- k[k >= 1L & k <= Km]; if (!length(k)) return(invisible(NULL))
    if (is.null(w)) w <- rep(1, length(k))
    f[[name]] <<- list(k = as.integer(k), w = w/sum(w))
  }
  for (K in c(2L, 3L, 4L, 5L, 6L, 8L, 10L, 14L, 20L, 28L, 40L))
    if (K <= Km) add(sprintf("dense_%02d", K), 1:K)
  add("dense_max", 1:Km)
  add("geo2", as.integer(unique(c(2^(0:6), Km)[c(2^(0:6), Km) <= Km])))
  add("geo3", as.integer(unique(c(3^(0:4), Km)[c(3^(0:4), Km) <= Km])))
  add("w_invsq", 1:Km, 1/(1:Km)^2)
  add("w_geom",  1:Km, 2^-(1:Km))
  add("rule12",  p_rule(n))
  f
}

## ------------------------- stage 1: the null side --------------------
## One profile at K_max serves every order and every family, because the
## statistics are nested partial sums of a single cross-moment array --
## the same property that makes the scan cost one pass in Section 2.6.
null_batch <- function(n, Km, B, batch) {
  set.seed(null_seed(n, batch))
  M <- matrix(0, R0, Km)
  for (r in seq_len(R0)) M[r, ] <- sp(runif(n), runif(n), Km, "poly", "first", B)
  M
}
pvals <- function(M, laws) {
  R <- nrow(M)
  sapply(seq_len(ncol(M)), function(k) 1 - findInterval(M[, k], laws[[k]])/length(laws[[k]]))
}
minw <- function(P, fam) do.call(pmin, as.data.frame(P[, fam$k, drop = FALSE] /
                                  rep(fam$w, each = nrow(P))))

build_null <- function(n) {
  fn <- null_file(n); got <- read_rds(fn)
  if (!is.null(got) && identical(got$R0, R0) && identical(got$seed, SEED)) return(got)
  Km <- K_max(n); B <- sb(n, Km, "poly")
  M1 <- null_batch(n, Km, B, 1L)
  M2 <- null_batch(n, Km, B, 2L)
  M3 <- null_batch(n, Km, B, 3L)
  laws <- lapply(seq_len(Km), function(k) sort(M1[, k]))
  qk   <- vapply(laws, function(l) l[ceiling((1 - ALPHA)*R0)], 0)   # single-order cut-offs
  P2 <- pvals(M2, laws); P3 <- pvals(M3, laws)
  fam <- families(Km, n)
  cal <- lapply(fam, function(f) {
    c05 <- as.numeric(quantile(minw(P2, f), ALPHA, type = 1))
    list(c_alpha = c05, size = mean(minw(P3, f) <= c05))
  })
  ## moment and correlation diagnostics: where the profile stops behaving
  ## like a chi-squared partial sum, which is the information boundary.
  cr  <- cor(M1)
  dev <- vapply(seq_len(Km), function(k) {
    if (k == Km) return(NA_real_)
    j <- seq_len(Km)[-k]
    max(abs(cr[k, j] - pmin(k, j)/pmax(k, j)))
  }, 0)
  mom <- data.frame(n = n, k = seq_len(Km), k2_over_n = seq_len(Km)^2/n,
                    mean_ratio = colMeans(M1)/seq_len(Km)^2,
                    var_ratio  = apply(M1, 2, var)/(2*seq_len(Km)^2),
                    corr_maxdev = dev,
                    size_single = colMeans(P3 <= ALPHA))
  out <- list(n = n, Km = Km, R0 = R0, seed = SEED, laws = laws, qk = qk,
              fam = fam, cal = cal, mom = mom)
  save_rds(out, fn)
  out
}

message("\n--- stage 1: null calibration ---")
NULLS <- list()
for (n in N_GRID) {
  t0 <- Sys.time(); NULLS[[as.character(n)]] <- build_null(n)
  message(sprintf("  n = %4d  K_max = %2d  %5.0f s", n, K_max(n),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
}

## ------------------------ stage 2: the power side --------------------
## Returns, for one cell, the rejection rate of every single order and of
## every family, all on the same R1 datasets.
run_cell <- function(cell) {
  fn <- cell_file(cell$n, cell$si, cell$nu)
  if (!is.null(read_rds(fn))) return(invisible(NULL))
  n <- cell$n; nl <- NULLS[[as.character(n)]]; Km <- nl$Km
  B <- sb(n, Km, "poly"); g <- gen[[cell$si]]
  set.seed(cell_seed(n, cell$si, cell$nu))
  M <- matrix(0, R1, Km)
  for (r in seq_len(R1)) {
    d <- g(n, cell$nu)
    M[r, ] <- sp(d$x, d$y, Km, "poly", "first", B)
  }
  pw <- vapply(seq_len(Km), function(k) mean(M[, k] > nl$qk[k]), 0)
  P  <- pvals(M, nl$laws)
  fp <- vapply(names(nl$fam), function(nm)
         mean(minw(P, nl$fam[[nm]]) <= nl$cal[[nm]]$c_alpha), 0)
  save_rds(list(n = n, si = cell$si, scheme = names(gen)[cell$si], nu = cell$nu,
                order_power = pw, family_power = fp), fn)
  invisible(NULL)
}

cells <- list()
for (n in N_GRID) for (si in seq_along(gen)) for (nu in NU)
  cells[[length(cells) + 1L]] <- list(n = n, si = si, nu = nu)
todo <- Filter(function(cl) is.null(read_rds(cell_file(cl$n, cl$si, cl$nu))), cells)
message(sprintf("\n--- stage 2: power ---\n  %d cells, %d to run, %d cores",
                length(cells), length(todo), NCORES))

## A cell costs about n K_max(n)^2, and with K_max = floor(sqrt(n)) that is
## of order n^2: an n = 1600 cell is some sixty times an n = 25 one.  Timing
## the first cell of the grid and multiplying by the number of cells, as an
## earlier version did, therefore understated the run by a factor of five.
## Time one cell at EVERY sample size instead -- they are cells that have to
## be computed anyway, so nothing is spent on the measurement -- and project
## the rest from the measured cost of its own n.  This is the lesson already
## paid for in 02_power_continuous.R.
if (length(todo)) {
  secs <- setNames(rep(NA_real_, length(N_GRID)), as.character(N_GRID))
  first <- vapply(N_GRID, function(n) {
    w <- which(vapply(todo, function(cl) cl$n == n, TRUE)); if (length(w)) w[1] else NA_integer_
  }, 0L)
  for (i in which(!is.na(first))) {
    t0 <- Sys.time(); run_cell(todo[[first[i]]])
    secs[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    message(sprintf("  n = %4d: one cell takes %6.1f s", N_GRID[i], secs[i]))
  }
  ## Any sample size already finished before this run gets the cost of the
  ## nearest measured one, scaled by n K_max(n)^2, so a resumed run still
  ## projects sensibly.
  wgt <- function(n) n*K_max(n)^2
  base <- which(!is.na(secs))[1]
  cost <- function(n) {
    i <- match(as.character(n), names(secs))
    if (!is.na(secs[i])) return(secs[i])
    secs[base]*wgt(n)/wgt(N_GRID[base])
  }
  todo <- todo[-stats::na.omit(first)]
  if (length(todo)) {
    tot <- sum(vapply(todo, function(cl) cost(cl$n), 0))
    message(sprintf("  %d cells left, approx %.0f min on %d cores (%.0f min of it at n = %d)",
                    length(todo), tot/60/NCORES, NCORES,
                    sum(vapply(todo, function(cl)
                      if (cl$n == max(N_GRID)) cost(cl$n) else 0, 0))/60/NCORES,
                    max(N_GRID)))
  }
  if (length(todo)) {
    cl <- makeCluster(NCORES)
    clusterExport(cl, c("gen", "sp", "sb", "NULLS", "R1", "ALPHA", "SEED", "K_max",
                        "OUTDIR", "cell_seed", "cell_file", "save_rds", "read_rds",
                        "pvals", "minw", "p_rule", "run_cell"), envir = environment())
    invisible(clusterEvalQ(cl, library(dependence)))
    invisible(parLapplyLB(cl, todo, run_cell))
    stopCluster(cl)
  }
}

## ------------------------------ assemble -----------------------------
files <- list.files(OUTDIR, CELL_PAT, full.names = TRUE)
rows <- list(); ord <- list()
for (f in files) {
  z <- read_rds(f); if (is.null(z)) { warning("unreadable: ", f); next }
  rows[[length(rows) + 1L]] <- data.frame(
    n = z$n, scheme = z$scheme, nu = z$nu,
    family = names(z$family_power), power = as.numeric(z$family_power),
    oracle = max(z$order_power), best_k = which.max(z$order_power),
    rule12 = z$order_power[p_rule(z$n)], stringsAsFactors = FALSE)
  ord[[length(ord) + 1L]] <- data.frame(
    n = z$n, scheme = z$scheme, nu = z$nu, k = seq_along(z$order_power),
    power = z$order_power, stringsAsFactors = FALSE)
}
pow <- do.call(rbind, rows); ords <- do.call(rbind, ord)
pow$regret <- pow$oracle - pow$power
write.csv(pow,  paste0(OUTSTEM, "_power.csv"),  row.names = FALSE)
write.csv(ords, paste0(OUTSTEM, "_orders.csv"), row.names = FALSE)

nullrows <- do.call(rbind, lapply(NULLS, function(z) do.call(rbind, lapply(
  names(z$fam), function(nm) data.frame(
    n = z$n, family = nm, orders = paste(z$fam[[nm]]$k, collapse = "+"),
    n_orders = length(z$fam[[nm]]$k), top_k = max(z$fam[[nm]]$k),
    c_alpha = z$cal[[nm]]$c_alpha,
    ## The rule is min_k p_k / w_k <= c, so the union bound gives exactly
    ## c = ALPHA whatever the weights; the gain over it is therefore
    ## c_alpha/ALPHA, and the number of orders divided by that gain is the
    ## effective multiplicity, which for uniform weights is ALPHA over the
    ## unweighted min-p threshold.  Reporting the raw c_alpha alone would
    ## compare weighted and unweighted families on different scales.
    union_gain = z$cal[[nm]]$c_alpha/ALPHA,
    eff_mult = length(z$fam[[nm]]$k)*ALPHA/z$cal[[nm]]$c_alpha,
    size = z$cal[[nm]]$size, stringsAsFactors = FALSE)))))
write.csv(nullrows, paste0(OUTSTEM, "_null.csv"), row.names = FALSE)
moments <- do.call(rbind, lapply(NULLS, `[[`, "mom"))
write.csv(moments, paste0(OUTSTEM, "_moments.csv"), row.names = FALSE)

## ------------------------------- summary -----------------------------
## Regret is meaningless where every test rejects or none does, so the
## minimax is taken over the INFORMATIVE cells only, those in which the
## oracle itself lies strictly inside the band.  Cells outside it are
## reported separately so that the filter cannot hide a failure.
INFORM <- pow$oracle > 0.15 & pow$oracle < 0.95
con <- file(paste0(OUTSTEM, "_summary.txt"), "w")
say <- function(...) { cat(..., "\n", sep = "", file = con); cat(..., "\n", sep = "") }

say("K(n) calibration, seed ", SEED, ", R0 = ", R0, " x 3 batches, R1 = ", R1)
say(sum(INFORM), " informative cells of ", nrow(pow),
    " (oracle power strictly between 0.15 and 0.95)\n")

say("attained size on the third null batch, worst family per n (nominal ", ALPHA, "):")
for (n in N_GRID) {
  s <- nullrows[nullrows$n == n, ]
  say(sprintf("  n = %4d   min %.4f   max %.4f   (%s)", n, min(s$size), max(s$size),
              s$family[which.max(abs(s$size - ALPHA))]))
}
if (max(abs(nullrows$size - ALPHA))/ALPHA > 0.15)
  say("  WARNING: a family is off nominal by more than 15 per cent of it.",
      " R0 must be large enough that 1/R0 is far below c_alpha; at R0 = ", R0,
      " the p-value resolution is ", signif(1/R0, 2), ".")

say("\neffective multiplicity (orders / union-bound gain) against log K, dense families:")
d <- nullrows[grepl("^dense_[0-9]", nullrows$family), ]
d$ratio <- d$eff_mult/log(d$top_k)
for (n in N_GRID) {
  s <- d[d$n == n, ]
  if (nrow(s)) say(sprintf("  n = %4d   ratio in [%.2f, %.2f] over K = %s",
                           n, min(s$ratio), max(s$ratio), paste(s$top_k, collapse = ",")))
}

say("\nworst-case regret against the oracle order, over informative cells:")
## Row order: the family list of the largest n, so dense families come in
## order of K and the ladders and weighted families follow, with anything
## defined only at a smaller n appended rather than dropped.
fams <- names(NULLS[[as.character(max(N_GRID))]]$fam)
fams <- c(fams, setdiff(unique(pow$family), fams))
tab <- matrix(NA_real_, length(fams), length(N_GRID), dimnames = list(fams, N_GRID))
for (fm in fams) for (n in N_GRID) {
  s <- pow[INFORM & pow$family == fm & pow$n == n, ]
  if (nrow(s)) tab[fm, as.character(n)] <- max(s$regret)
}
say(paste(capture.output(print(round(tab, 3))), collapse = "\n"))

say("\nmean regret over informative cells:")
tab2 <- tab
for (fm in fams) for (n in N_GRID) {
  s <- pow[INFORM & pow$family == fm & pow$n == n, ]
  if (nrow(s)) tab2[fm, as.character(n)] <- mean(s$regret)
}
say(paste(capture.output(print(round(tab2, 3))), collapse = "\n"))

say("\nthe flat region: dense families within 0.01 of the best dense family,")
say("which is the interval of K over which the choice does not matter.")
for (n in N_GRID) {
  s <- tab[grepl("^dense_[0-9]", rownames(tab)), as.character(n)]
  s <- s[!is.na(s)]
  if (!length(s)) next
  flat <- names(s)[s <= min(s) + 0.01]
  say(sprintf("  n = %4d   best %s (%.3f)   flat over K = %s", n,
              names(s)[which.min(s)], min(s),
              paste(sub("dense_0?", "", flat), collapse = ", ")))
}

say("\nbest K per n, and the rules of the family of (12) that reproduce it:")
bestK <- vapply(N_GRID, function(n) {
  s <- tab[grepl("^dense_[0-9]", rownames(tab)), as.character(n)]; s <- s[!is.na(s)]
  if (!length(s)) return(NA_real_)
  as.numeric(sub("dense_0?", "", names(s)[which.min(s)]))
}, 0)
say("  n    ", paste(sprintf("%6d", N_GRID), collapse = ""))
say("  K*   ", paste(sprintf("%6.0f", bestK), collapse = ""))
ok <- !is.na(bestK)
if (sum(ok) >= 3) {
  grid <- expand.grid(r = seq(0.20, 0.55, by = 0.01), c = -2:2)
  grid$maxdev <- apply(grid, 1, function(g)
    max(abs(pmax(1, floor(N_GRID[ok]^g[["r"]]) + g[["c"]]) - bestK[ok])))
  grid <- grid[order(grid$maxdev), ]
  say("  closest rules K = max(1, floor(n^r) + c):")
  for (i in 1:5)
    say(sprintf("    r = %.2f, c = %+d   max deviation %.0f order(s)",
                grid$r[i], grid$c[i], grid$maxdev[i]))
}

say("\nthe oracle single order over informative cells: this is the direct")
say("evidence on the rate, since K must grow at least as fast as the order")
say("the best alternatives actually need.")
ok2 <- pow$family == fams[1] & INFORM          # one row per cell
say("  n      median k*   90th pct   max k*   K_max of grid")
for (n in N_GRID) {
  s <- pow[ok2 & pow$n == n, ]
  if (!nrow(s)) next
  say(sprintf("  %4d   %9.0f   %8.0f   %6.0f   %6d", n, median(s$best_k),
              quantile(s$best_k, 0.9, type = 1), max(s$best_k), K_max(n)))
}
say("  by scheme, median best order:")
for (sc in names(gen)) {
  s <- pow[ok2 & pow$scheme == sc, ]
  if (!nrow(s)) next
  m <- vapply(N_GRID, function(n) {
    z <- s$best_k[s$n == n]; if (length(z)) median(z) else NA_real_ }, 0)
  say(sprintf("    %-10s %s", sc, paste(sprintf("%5.0f", m), collapse = "")))
}

say("\nwhere the profile stops behaving like a chi-squared partial sum:")
for (n in N_GRID) {
  m <- moments[moments$n == n, ]
  k10 <- m$k[m$var_ratio > 1.10]; k25 <- m$k[m$var_ratio > 1.25]
  say(sprintf("  n = %4d   var/(2k^2) passes 1.10 at k = %s (k^2/n = %s), 1.25 at k = %s",
              n,
              if (length(k10)) k10[1] else "never",
              if (length(k10)) sprintf("%.2f", k10[1]^2/n) else "-",
              if (length(k25)) k25[1] else "never"))
}
close(con)
message("\nWritten: ", OUTSTEM, "_{power,orders,null,moments}.csv and _summary.txt")
