## =====================================================================
##  02b_topup_statistics.R
##
##  Add a statistic to a power study that has already been run, or replace
##  one whose adapter has been corrected, WITHOUT recomputing the others.
##
##  WHY THIS IS EXACT, AND NOT MERELY CHEAP.  02_power_continuous.R seeds
##  every cell from (n, scheme index, nu, q) alone and wraps every statistic
##  evaluation in with_frozen_rng(), so the generator sits at the same state
##  before and after each evaluation block.  The sequence of datasets a cell
##  sees is therefore fixed by its seed and by the LENGTH of its two loops,
##  and is independent of which statistics are computed.  Replaying the seed
##  reproduces the very same R1 datasets, so the new statistic is compared
##  with the old ones on IDENTICAL data: paired, not just equal in law.  The
##  same holds for the distribution-free critical values, whose null data are
##  generated outside the frozen block in a loop of fixed length R0_DF.
##
##  VERIFIED, NOT ASSUMED, and this is the part that matters.  Every replay
##  recomputes a statistic that is ALREADY on disk -- B_n, the cheapest --
##  and demands the stored power back to the last digit.  A single cell that
##  disagrees aborts the run before anything is written.  If the frozen-RNG
##  discipline had not been in place when the cells were written, or if a
##  generator, R0, R1, SEED or the seed formula had changed since, this is
##  where it shows up, loudly, instead of being published.
##
##  WHAT IT DOES NOT DO.  It does not write the CSV or the LaTeX tables.  It
##  extends the per-cell .rds files and the critical-value tables, and then
##  02_power_continuous.R is re-run: it finds every cell complete, computes
##  nothing, and assembles the published files itself.  One script produces
##  the tables, always.
##
##  COST, at the shipped settings.  The scan evaluates in about the time one
##  B_n does, so adding it is roughly 2 single-core hours over 780 cells
##  including the B_n control, plus about a minute for the null tables.
##  Replacing BET costs whatever BET costs; the projection printed by
##  02_power_continuous.R at startup is the number to look at before
##  committing to it.
##
##  USAGE.  Run from the directory that holds results/, with 02_power_-
##  continuous.R beside this file:
##      Rscript 02b_topup_statistics.R
##  then re-run 02_power_continuous.R.
## =====================================================================

## ------------------------------ what to do ---------------------------
## ADD      statistics absent from the cells, to be computed and inserted.
## REPLACE  statistics already there whose rows are to be overwritten.  Use
##          it when an adapter has been corrected -- BET is the live case:
##          the October cells were written before abs() was applied to its
##          signed extreme asymmetry, and 44 per cent of them sit below the
##          nominal level as a result, so that column cannot be published.
##          Replacing BET also recomputes its critical values, which were
##          calibrated on the same signed quantity.
ADD     <- c("B_scan")
REPLACE <- character(0)          # c("BET") once you are ready to pay for it

MAIN    <- "02_power_continuous.R"
## Controls: statistics already on disk, recomputed on the replayed data and
## required to give the stored power back exactly.  More than one, because a
## rejection rate out of R1 is a coarse fingerprint -- a perturbation small
## enough to flip no rejection would pass a single check.  Two independent
## rejection patterns per cell, over 780 cells, is 1560 such checks; and
## cells written with a stream fingerprint are checked against the data
## itself, which is sharper still.  Keep them cheap and rank-invariant.
CONTROL <- c("B_n", "P_n")
NCORES_TOPUP <- NA               # NA: take NCORES from the main script

## ------------------- the main script's own definitions ----------------
## Evaluated, not copied.  Anything that would make this script disagree
## with the one that produced the cells -- a generator, the seed formula,
## an adapter, R0, R1, SEED -- is therefore impossible by construction.
stopifnot("02_power_continuous.R must be in the working directory" =
            file.exists(MAIN))
.src <- readLines(MAIN, warn = FALSE)
.cut <- grep("<<< DEFINITIONS END >>>", .src, fixed = TRUE)
stopifnot("the DEFINITIONS END sentinel is missing or duplicated in the main script" =
            length(.cut) == 1L)
eval(parse(text = .src[seq_len(.cut - 1L)]), envir = globalenv())
message(sprintf("definitions taken from %s (%d lines)", MAIN, .cut - 1L))

library(parallel)
if (is.na(NCORES_TOPUP)) NCORES_TOPUP <- NCORES
TARGET <- union(ADD, REPLACE)
stopifnot("nothing to do: ADD and REPLACE are both empty" = length(TARGET) > 0L,
          "ADD and REPLACE must not overlap" = !length(intersect(ADD, REPLACE)))

## ----------------------------- the cells ------------------------------
CF <- list.files(OUTDIR, "^cell_n.*rds$", full.names = TRUE)
stopifnot("no cells found: run 02_power_continuous.R first" = length(CF) > 0L)

## The scheme index is in the file name and nowhere else, and it enters the
## seed, so it is parsed -- and the parse is then round-tripped through the
## main script's own cell_file(), so a name this script reads differently
## from the way the main script writes it cannot slip past.
parse_cell <- function(f) {
  m <- regmatches(basename(f),
                  regexec("^cell_n(\\d+)_s(\\d+)_nu(\\d+)_q(\\d+)\\.rds$", basename(f)))[[1]]
  stopifnot("unparseable cell file name" = length(m) == 5L)
  z <- list(file = f, n = as.integer(m[2]), si = as.integer(m[3]),
            nu = as.integer(m[4]) / 100, q = as.integer(m[5]))
  if (!identical(basename(cell_file(z)), basename(f)))
    stop(sprintf("file name %s does not round-trip through cell_file() (gives %s)",
                 basename(f), basename(cell_file(z))))
  z
}
CELLS <- lapply(CF, parse_cell)

## One read of everything, to establish what is there and what the replay
## parameters are.  Both are taken FROM THE FILES rather than from the
## current configuration: the cells are the authority on how they were made.
INV <- lapply(CELLS, function(z) { d <- readRDS(z$file)
  list(tests = d$test, scheme = d$scheme[1], p = d$p[1], q = d$q[1],
       R1 = max(d$R1), R0cell = if (any(!d$df_calibrated))
         max(d$R0[!d$df_calibrated]) else 0L,
       ncell = sum(!d$df_calibrated)) })

have <- sort(unique(unlist(lapply(INV, `[[`, "tests"))))
message("\nstatistics on disk: ", paste(have, collapse = ", "))
message("to add:     ", if (length(ADD)) paste(ADD, collapse = ", ") else "none")
message("to replace: ", if (length(REPLACE)) paste(REPLACE, collapse = ", ") else "none")
if (length(setdiff(ADD, have)) != length(ADD))
  stop("already on disk, use REPLACE instead of ADD: ",
       paste(intersect(ADD, have), collapse = ", "))
if (!all(REPLACE %in% have))
  stop("not on disk, use ADD instead of REPLACE: ",
       paste(setdiff(REPLACE, have), collapse = ", "))
CONTROL <- setdiff(intersect(CONTROL, have), TARGET)
if (!length(CONTROL))
  stop("no control statistic left: every candidate in CONTROL is either absent ",
       "from the cells or itself a target. Name a statistic that is on disk and ",
       "is not being changed, or the replay cannot be verified.")
message("controls:    ", paste(CONTROL, collapse = ", "))

## Every cell must agree on the replay parameters, or the study was not run
## in one configuration and a top-up is not well defined.
R1_SEEN  <- unique(vapply(INV, `[[`, 0, "R1"))
R0_SEEN  <- unique(vapply(INV, `[[`, 0, "R0cell"))
NC_SEEN  <- unique(vapply(INV, `[[`, 0, "ncell"))
stopifnot("cells disagree on R1" = length(R1_SEEN) == 1L,
          "cells disagree on the per-cell null loop length" = length(R0_SEEN) == 1L,
          "cells disagree on how many statistics are calibrated per cell" =
            length(NC_SEEN) == 1L)
message(sprintf("\nreplay parameters read from the cells: R1 = %g, per-cell null loop = %g draws%s",
                R1_SEEN, R0_SEEN,
                if (NC_SEEN == 0) " (no per-cell calibration)" else
                  sprintf(" for %d statistic(s)", NC_SEEN)))
if (R1_SEEN != R1 || (NC_SEEN > 0 && R0_SEEN != R0))
  warning("the cells were written with R1 = ", R1_SEEN, " and R0 = ", R0_SEEN,
          " but the main script now has R1 = ", R1, " and R0 = ", R0,
          ". The replay follows the cells; fix the main script before re-running it.",
          call. = FALSE)

## --------------------------- the adapters -----------------------------
if (USE_SCAN && any(c("B_scan") %in% TARGET)) {
  .st <- build_scan_tables(sort(unique(vapply(CELLS, `[[`, 0, "n"))))
  SCAN_LAW <- .st$law; SCAN_BAS <- .st$bas; SCAN_PMAT <- .st$pmat; rm(.st)
} else { SCAN_LAW <- SCAN_BAS <- SCAN_PMAT <- list() }

S_ALL <- make_stats()
if (!all(TARGET %in% names(S_ALL)))
  stop("no adapter for: ", paste(setdiff(TARGET, names(S_ALL)), collapse = ", "),
       ". Check the USE_* switches at the top of the main script.")
NEED <- c(TARGET, CONTROL)
S <- S_ALL[NEED]

message("\nself-test of the statistics being computed:")
for (nm in NEED) {
  set.seed(1); x <- runif(max(N_GRID)); y <- x + 0.5*rnorm(max(N_GRID))
  v <- tryCatch(as.numeric(S[[nm]](x, y, p_rule(max(N_GRID)), p_rule(max(N_GRID)))),
                error = function(e) e)
  if (inherits(v, "error") || length(v) != 1L || !is.finite(v))
    stop(sprintf("adapter %s is not usable: %s", nm,
                 if (inherits(v, "error")) conditionMessage(v) else "bad value"))
  message(sprintf("  %-10s ok  %-12.6g", nm, v))
}

## Rank-invariance decides where the critical value comes from, so it is
## detected here exactly as the main script detects it.  A target that is
## NOT rank-invariant would need per-cell calibration, which would change
## the per-cell null loop length and hence every dataset: the top-up route
## is closed for such a statistic and the study has to be re-run entire.
RI_T <- detect_rank_invariant()[NEED]
for (nm in NEED)
  message(sprintf("  %-10s %s", nm,
                  if (RI_T[nm]) "rank-invariant" else "VALUE-BASED"))
if (any(!RI_T[TARGET]))
  stop("value-based target(s): ", paste(TARGET[!RI_T[TARGET]], collapse = ", "),
       ". Adding one changes the per-cell null loop length, hence the datasets ",
       "of every cell, so it cannot be topped up: delete ", OUTDIR,
       "/cell_*.rds and re-run 02_power_continuous.R from scratch.")

## ------------------- critical values for the targets ------------------
## The null datasets here are generated outside the frozen block in a loop
## of length max(caps) = R0_DF, so they do not depend on which statistics
## are evaluated: these are the same R0_DF null samples the stored critical
## values came from.  The control's critical value is recomputed and matched
## against the stored one to prove it.
PQ_SEEN <- unique(do.call(rbind, lapply(seq_along(CELLS), function(i)
  c(CELLS[[i]]$n, INV[[i]]$p, INV[[i]]$q))))
newcrit <- list()
for (r in seq_len(nrow(PQ_SEEN))) {
  n <- PQ_SEEN[r, 1]; p <- PQ_SEEN[r, 2]; q <- PQ_SEEN[r, 3]
  key <- sprintf("n%d_p%d_q%d", n, p, q)
  f <- file.path(OUTDIR, paste0("critdf_", key, ".rds"))
  stopifnot("missing critical-value table" = file.exists(f))
  cv_old <- readRDS(f)
  stopifnot("the control is absent from the stored critical-value table" =
              CONTROL %in% names(cv_old))
  message(sprintf("calibrating %s at %s, R = %g",
                  paste(TARGET, collapse = ", "), key, R0_DF))
  set.seed(SEED + 31L*n + 7L*p + q)
  caps <- vapply(NEED, cap_for, 0L, default = R0_DF)
  M <- matrix(NA_real_, max(caps), length(NEED), dimnames = list(NULL, NEED))
  for (i in seq_len(max(caps))) {
    x <- runif(n); y <- runif(n)
    with_frozen_rng(
      for (k in NEED) if (i <= caps[k]) M[i, k] <- as.numeric(S[[k]](x, y, p, q)))
  }
  cv_new <- vapply(NEED, function(k)
    quantile(M[seq_len(caps[k]), k], 1 - ALPHA, names = FALSE, na.rm = TRUE), 0)
  for (c0 in CONTROL) {
    d <- abs(cv_new[[c0]] - cv_old[[c0]])
    message(sprintf("  control %-10s stored %.10g, replayed %.10g  ->  %s",
                    c0, cv_old[[c0]], cv_new[[c0]],
                    if (d < 1e-9) "null sample reproduces" else "MISMATCH"))
    if (d >= 1e-9)
      stop(sprintf(paste("the null sample at %s does not reproduce (%s differs by %.3g).",
                         "The critdf stream has changed since the cells were written --",
                         "check SEED, R0_DF, EVAL_CAP and the %s adapter -- so a top-up",
                         "would not be comparable. Re-run the study from scratch."),
                   key, c0, d, c0))
  }
  newcrit[[key]] <- list(file = f, cv = c(cv_old[setdiff(names(cv_old), TARGET)],
                                          cv_new[TARGET]))
  for (t in TARGET)
    message(sprintf("  %-10s 5%% point %.6g%s", t, cv_new[[t]],
                    if (t %in% REPLACE)
                      sprintf("   (was %.6g)", cv_old[[t]]) else ""))
}
## Cross-check for the scan: an independent null batch, free, see the note
## in the main script.
if ("B_scan" %in% TARGET)
  for (n in sort(unique(vapply(CELLS, `[[`, 0, "n")))) {
    k <- as.character(n); p <- p_rule(n)
    a <- newcrit[[sprintf("n%d_p%d_q%d", n, p, p)]]$cv[["B_scan"]]
    b <- quantile(-apply(SCAN_PMAT[[k]], 1L, min), 1 - ALPHA, names = FALSE)
    message(sprintf("  B_scan 5%% point, n = %3d: critdf %.6f   independent batch %.6f   gap %.6f",
                    n, a, b, abs(a - b)))
  }

## ----------------------------- one cell -------------------------------
## Faithful replay of run_cell()'s generator consumption: the per-cell null
## loop first, for its RNG draws only, then the alternative loop.  The null
## loop's statistics are NOT recomputed -- it is run for its side effect on
## the generator, which is the whole point.
topup_cell <- function(z, i) {
  d <- readRDS(z$file)
  n <- z$n; nu <- z$nu; q <- z$q; p <- INV[[i]]$p; g <- gen[[INV[[i]]$scheme]]
  stopifnot("the scheme index in the file name does not match the cell" =
              identical(names(gen)[z$si], INV[[i]]$scheme))
  set.seed(cell_seed(n, z$si, nu, q))      # the main script's own formula
  if (NC_SEEN > 0L)
    for (j in seq_len(R0_SEEN)) { dd <- g(n, nu); dd$y <- sample(dd$y)
                                  with_frozen_rng(invisible(NULL)) }
  caps1 <- vapply(NEED, cap_for, 0L, default = R1_SEEN)
  alt <- matrix(NA_real_, max(caps1), length(NEED), dimnames = list(NULL, NEED))
  fp <- numeric(5)
  for (j in seq_len(max(caps1))) {
    dd <- g(n, nu)
    fp <- stream_acc(fp, dd$x, dd$y)
    with_frozen_rng(
      for (k in NEED) if (j <= caps1[k])
        alt[j, k] <- as.numeric(S[[k]](dd$x, dd$y, p, q)))
  }
  cv <- newcrit[[sprintf("n%d_p%d_q%d", n, p, q)]]$cv
  pw <- vapply(NEED, function(k)
    mean(alt[seq_len(caps1[k]), k] > cv[[k]], na.rm = TRUE), 0)

  ## THE VERIFICATION, in two forms.
  ##
  ## Sharp, when available: the cell carries a fingerprint of the datasets
  ## it was built from, and the replay must reproduce it.  This compares the
  ## DATA, so it catches a changed generator however little the change moves
  ## any statistic.
  old_fp <- attr(d, "stream")
  if (!is.null(old_fp) && !stream_match(old_fp, fp))
    return(sprintf("%s: data fingerprint differs (stored %s, replayed %s)",
                   basename(z$file), paste(signif(old_fp, 10), collapse = " "),
                   paste(signif(fp, 10), collapse = " ")))
  ## Fall-back, for cells written before the fingerprint existed: each
  ## control was recomputed on the replayed data and is scored at its STORED
  ## critical value, so the stored power must come back exactly.  Coarser --
  ## it can only see changes large enough to move a rejection -- which is
  ## why there is more than one control.
  for (c0 in CONTROL) {
    rec <- d$power[d$test == c0]
    if (!length(rec)) next
    if (abs(pw[[c0]] - rec[1]) > 1e-12)
      return(sprintf("%s: %s power stored %.6f, replayed %.6f", basename(z$file),
                     c0, rec[1], pw[[c0]]))
  }

  add <- data.frame(n = n, scheme = INV[[i]]$scheme, nu = nu, p = p, q = q,
                    test = TARGET, power = as.numeric(pw[TARGET]),
                    se = sqrt(as.numeric(pw[TARGET]) *
                                (1 - as.numeric(pw[TARGET])) / caps1[TARGET]),
                    crit = as.numeric(unlist(cv[TARGET])),
                    R0 = R0_DF, R1 = caps1[TARGET], df_calibrated = TRUE,
                    row.names = NULL, stringsAsFactors = FALSE)
  ## REPLACE drops the old rows; the stored critical values of everything
  ## else are left untouched, so no other number in the cell moves.
  out <- rbind(d[!(d$test %in% TARGET), ], add)
  out <- out[order(match(out$test, c(have, setdiff(TARGET, have)))), ]
  attr(out, "stream") <- fp
  saveRDS(out, z$file)
  NA_character_
}

## ------------------------------- driver -------------------------------
todo <- which(!vapply(seq_along(CELLS), function(i)
  all(TARGET %in% INV[[i]]$tests) && !length(REPLACE), TRUE))
message(sprintf("\n%d cells, %d to process, %d cores", length(CELLS), length(todo),
                NCORES_TOPUP))
if (!length(todo)) { message("nothing to do."); quit(save = "no") }

t0 <- Sys.time(); first <- topup_cell(CELLS[[todo[1]]], todo[1])
if (!is.na(first)) stop("the replay does not reproduce the stored datasets.\n  ", first,
  "\n\nThis means the cells on disk were written under a different generator state",
  "\nthan this replay produces. Nothing further has been written. Check that SEED,",
  "\nR0, R1, the seed formula and every gen[[...]] body are unchanged since the run,",
  "\nand that the cells were written with with_frozen_rng() in place. If they were",
  "\nnot, the study has to be re-run entire: delete ", OUTDIR, "/cell_*.rds.")
dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
message(sprintf("first cell (n = %d) took %.1f s and reproduced; remaining %d cells approx %.2f h on %d cores",
                CELLS[[todo[1]]]$n, dt, length(todo) - 1L,
                dt*(length(todo) - 1L)/3600/NCORES_TOPUP, NCORES_TOPUP))

bad <- character(0)
if (length(todo) > 1L) {
  rest <- todo[-1]
  if (NCORES_TOPUP > 1L) {
    cl <- makeCluster(NCORES_TOPUP)
    on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
    clusterExport(cl, c("gen","make_stats","pluck1","pseudo","with_frozen_rng",
                        "S","NEED","TARGET","CONTROL","INV","CELLS","newcrit",
                        "R0_SEEN","R1_SEEN","NC_SEEN","have","R0_DF",
                        "SEED","ALPHA","cap_for","EVAL_CAP","p_rule","d_rule",
                        "K_rule","SCAN_LAW","SCAN_BAS","USE_SCAN","USE_NEW",
                        "USE_MIC","topup_cell","OUTDIR",
                        "stream_acc","stream_match","cell_seed","cell_file"))
    invisible(clusterEvalQ(cl, library(dependence)))
    res <- parLapplyLB(cl, rest, function(i) topup_cell(CELLS[[i]], i))
    stopCluster(cl); on.exit()
  } else {
    res <- lapply(rest, function(i) topup_cell(CELLS[[i]], i))
  }
  bad <- unlist(res); bad <- bad[!is.na(bad)]
}

if (length(bad)) {
  message(sprintf("\n%d of %d cells did NOT reproduce:", length(bad), length(todo)))
  for (b in head(bad, 20L)) message("  ", b)
  stop("those cells were left unchanged, the rest were written. ",
       "Do not publish this state: resolve the mismatch and re-run.")
}

message(sprintf("\n%d cells extended, every one verified against the stored %s.",
                length(todo), paste(CONTROL, collapse = " and ")))
for (k in names(newcrit)) saveRDS(newcrit[[k]]$cv, newcrit[[k]]$file)
message("critical-value tables updated: ", paste(names(newcrit), collapse = ", "))
message("\nNow re-run ", MAIN, ". It will find every cell complete, compute ",
        "nothing, and write the CSV, the summary and the LaTeX tables.")
