## =====================================================================
##  pilot_multibasis.R -- EXPLORATORY PILOT, not part of the paper's grid
##
##  Question: does combining basis families that look in different
##  directions improve the power of the rank canonical-correlation test,
##  and at what cost where a single family is already right?
##
##  PRE-REGISTERED DESIGN.  Everything in this header was fixed before any
##  power was computed, and none of it was changed afterwards.
##
##  1. The alternatives are deliberately DISJOINT from the thirteen schemes
##     of the comparison study, so that whatever is learned here can be
##     taken to that study as an out-of-sample check rather than having
##     been tuned on it.  Five, chosen to span frequency and locality:
##       gauss_copula  Gaussian copula, rho = 0.05..0.50   monotone, lowest
##                     frequency: the case where one family is already right
##                     and any combination can only cost
##       chirp         Y = sin(8 pi X^2) + noise           frequency rising
##                     from 0 to 8 cycles across the range of X
##       square        Y = sign(sin(6 pi X)) + noise       discontinuous:
##                     3 cycles plus odd harmonics beyond the octaves used
##       bump          Y = exp(-((X - 0.7)/0.05)^2) + noise   localised in X
##       trig_copula   copula density 1 + theta sin(8 pi u) sin(8 pi v),
##                     theta = 0.1..1   joint oscillation in both variables
##     The noise of the three functional schemes has sd nu * sd(f), with
##     nu = 0.2..2.0, so that nu = 1 is unit signal-to-noise in each.
##
##  2. The candidate tests, all built from P_n on grid-orthonormal bases of
##     the ranks (B_n behaves identically to within 0.003 in mean power):
##       fixed        P_n, polynomials, p = q from rule (12)
##       scan         the order scan of Section 2.6, K = p + 2
##       oct_fun      min-p over the fixed statistic and FUNCTIONAL octave
##                    pairs: cosine block x polynomials of order 2, both
##                    orientations, for octaves [4,7] and [8,15]
##       oct_sym      min-p over the fixed statistic and SYMMETRIC octave
##                    pairs: [4,7] x [4,7] and [8,15] x [8,15]
##       oct_both     min-p over both kinds of octave pair
##       oct_fun_sum  the SUM of the oct_fun statistics, not their minimum
##       haar_fun     as oct_fun with Haar levels 2 and 3 in place of the
##                    cosine octaves
##       scan_oct     the scan with the functional octave pairs added: the
##                    generalisation of Section 2.6 to several families
##
##     Why these, on principle and not on the schemes:
##     - Octave blocks [2^j, 2^(j+1) - 1] are equal steps in log-frequency,
##       the spacing the OU structure of the null profile singles out.
##       They start at 4 because cosines 1..3 span almost the same space as
##       the polynomials (overlap 2.74 of 3 at n = 128): adding them would
##       cost multiplicity and buy nothing.
##     - The functional pairing (oscillating basis on one variable, order-2
##       polynomials on the other) follows from Proposition 1 applied to
##       the regression model Y = f(X) + e: the local direction of
##       departure is f(u) psi(v), with psi the noise score, which loads
##       almost entirely on the first polynomial; order 2 adds scale.
##       Both orientations are included so the test stays symmetric in X, Y.
##     - The symmetric pairing covers departures in which both variables
##       oscillate together, which the functional pairing does not.
##
##  3. Calibration, as everywhere in this project: three independent null
##     batches -- per-component laws, combination thresholds, and a size
##     check on data used for neither.  All tests are size-adjusted at 5%.
##
##  4. Every variant is evaluated on the same datasets, so comparisons
##     between variants are paired.
## =====================================================================

if (dir.exists("/home/claude/pkgtest/lib")) .libPaths(c("/home/claude/pkgtest/lib", .libPaths()))
library(dependence)
library(parallel)

SEED  <- 20261016L
NS    <- c(100L, 200L)
R0    <- 20000L
R1    <- 2000L
ALPHA <- 0.05
## mclapply forks, which Windows cannot do: run serially there (a few minutes).
NCORE <- if (.Platform$OS.type == "windows") 1L else max(1L, detectCores())
p_rule <- function(n) max(1L, as.integer(floor(n^0.31)) - 1L)

## ------------------------------ bases --------------------------------
ugrid  <- function(n) (seq_len(n) - 0.5)/n
orthon <- function(M, n) {                       # centre, orthonormalise, scale
  M <- sweep(as.matrix(M), 2, colMeans(M))
  qr.Q(qr(M)) * sqrt(n)
}
make_bases <- function(n) {
  p <- p_rule(n); K <- p + 2L
  PK <- dependence:::scan_basis(n, K, "poly")[[K]]     # nested, crossprod = n I
  u  <- ugrid(n)
  haar <- function(l) orthon(sapply(0:(2^l - 1), function(m) {
    a <- m/2^l; b <- (m + 0.5)/2^l; e <- (m + 1)/2^l
    ifelse(u >= a & u < b, 1, ifelse(u >= b & u < e, -1, 0)) }), n)
  B <- list(n = n, p = p, K = K, PK = PK,
            C47  = orthon(sapply(4:7,  function(j) cos(pi*j*u)), n),
            C815 = orthon(sapply(8:15, function(j) cos(pi*j*u)), n),
            H2 = haar(2), H3 = haar(3))
  for (nm in c("PK", "C47", "C815", "H2", "H3")) {
    G <- crossprod(B[[nm]])/n
    stopifnot(max(abs(G - diag(ncol(G)))) < 1e-8, max(abs(colMeans(B[[nm]]))) < 1e-8)
  }
  B
}

## ---------------------- component statistics -------------------------
## One call per dataset returns every component that any variant uses, so
## all variants see the same data.  st(A, B) = n || A'B / n ||_F^2, the
## Pillai statistic for spans A and B (exact, since A'A = B'B = n I).
comp_names <- function(K) c(paste0("prof", seq_len(K)),
  "c47_p2", "p2_c47", "c815_p2", "p2_c815", "c47_c47", "c815_c815",
  "h2_p2", "p2_h2", "h3_p2", "p2_h3")
components <- function(rx, ry, B) {
  n <- B$n; K <- B$K
  st <- function(A, Bm) sum(crossprod(A, Bm)^2)/n
  PKx <- B$PK[rx, , drop = FALSE]; PKy <- B$PK[ry, , drop = FALSE]
  S2  <- (crossprod(PKx, PKy)/n)^2
  prof <- vapply(seq_len(K), function(k) n*sum(S2[seq_len(k), seq_len(k)]), 0)
  P2x <- PKx[, 1:2]; P2y <- PKy[, 1:2]
  Cx4 <- B$C47[rx, ];  Cy4 <- B$C47[ry, ]
  Cx8 <- B$C815[rx, ]; Cy8 <- B$C815[ry, ]
  Hx2 <- B$H2[rx, ];   Hy2 <- B$H2[ry, ]
  Hx3 <- B$H3[rx, ];   Hy3 <- B$H3[ry, ]
  c(prof, st(Cx4, P2y), st(P2x, Cy4), st(Cx8, P2y), st(P2x, Cy8),
    st(Cx4, Cy4), st(Cx8, Cy8), st(Hx2, P2y), st(P2x, Hy2), st(Hx3, P2y), st(P2x, Hy3))
}

variants <- function(p, K) {
  pp <- paste0("prof", p)
  fun <- c("c47_p2", "p2_c47", "c815_p2", "p2_c815")
  list(fixed       = list(type = "single", comp = pp),
       scan        = list(type = "minp",   comp = paste0("prof", seq_len(K))),
       oct_fun     = list(type = "minp",   comp = c(pp, fun)),
       oct_sym     = list(type = "minp",   comp = c(pp, "c47_c47", "c815_c815")),
       oct_both    = list(type = "minp",   comp = c(pp, fun, "c47_c47", "c815_c815")),
       oct_fun_sum = list(type = "sum",    comp = c(pp, fun)),
       haar_fun    = list(type = "minp",   comp = c(pp, "h2_p2", "p2_h2", "h3_p2", "p2_h3")),
       scan_oct    = list(type = "minp",   comp = c(paste0("prof", seq_len(K)), fun)))
}

## ----------------------------- the null ------------------------------
null_batch <- function(B, batch, R = R0) {
  set.seed(SEED + 1000L*B$n + 7L*batch)
  M <- t(vapply(seq_len(R), function(i)
         components(sample.int(B$n), sample.int(B$n), B), numeric(length(comp_names(B$K)))))
  colnames(M) <- comp_names(B$K); M
}
pvals <- function(M, laws) {
  P <- sapply(colnames(M), function(nm)
         1 - findInterval(M[, nm], laws[[nm]])/length(laws[[nm]]))
  matrix(P, nrow = nrow(M), dimnames = list(NULL, colnames(M)))
}
calibrate <- function(B) {
  M1 <- null_batch(B, 1L); M2 <- null_batch(B, 2L); M3 <- null_batch(B, 3L)
  laws <- lapply(colnames(M1), function(nm) sort(M1[, nm])); names(laws) <- colnames(M1)
  P2 <- pvals(M2, laws); P3 <- pvals(M3, laws)
  V <- variants(B$p, B$K)
  cal <- lapply(names(V), function(v) {
    s <- V[[v]]
    if (s$type == "single") {
      cut <- laws[[s$comp]][ceiling((1 - ALPHA)*R0)]
      list(cut = cut, size = mean(M3[, s$comp] > cut), meff = 1)
    } else if (s$type == "sum") {
      sm1 <- sort(rowSums(M1[, s$comp])); cut <- sm1[ceiling((1 - ALPHA)*R0)]
      list(cut = cut, size = mean(rowSums(M3[, s$comp]) > cut), meff = NA)
    } else {
      mp2 <- do.call(pmin, as.data.frame(P2[, s$comp, drop = FALSE]))
      cut <- as.numeric(quantile(mp2, ALPHA, type = 1))
      mp3 <- do.call(pmin, as.data.frame(P3[, s$comp, drop = FALSE]))
      list(cut = cut, size = mean(mp3 <= cut), meff = ALPHA/cut)
    }
  })
  names(cal) <- names(V)
  list(laws = laws, cal = cal, V = V, nullcor = cor(M1))
}

## -------------------------- the alternatives -------------------------
f_chirp  <- function(x) sin(8*pi*x^2)
f_square <- function(x) sign(sin(6*pi*x))
f_bump   <- function(x) exp(-((x - 0.7)/0.05)^2)
xx <- (seq_len(1e6) - 0.5)/1e6
SD <- c(chirp = sd(f_chirp(xx)), square = sd(f_square(xx)), bump = sd(f_bump(xx)))
gen <- list(
  gauss_copula = function(n, s) { z <- rnorm(n); list(x = z, y = s*z + sqrt(1 - s^2)*rnorm(n)) },
  chirp  = function(n, s) { x <- runif(n); list(x = x, y = f_chirp(x)  + s*SD[["chirp"]]*rnorm(n)) },
  square = function(n, s) { x <- runif(n); list(x = x, y = f_square(x) + s*SD[["square"]]*rnorm(n)) },
  bump   = function(n, s) { x <- runif(n); list(x = x, y = f_bump(x)   + s*SD[["bump"]]*rnorm(n)) },
  trig_copula = function(n, s) {
    out <- matrix(numeric(0), 0, 2)
    while (nrow(out) < n) {
      m <- 3L*n; u <- runif(m); v <- runif(m)
      keep <- runif(m) < (1 + s*sin(8*pi*u)*sin(8*pi*v))/(1 + s)
      out <- rbind(out, cbind(u, v)[keep, , drop = FALSE])
    }
    list(x = out[seq_len(n), 1], y = out[seq_len(n), 2])
  })
STRENGTH <- list(gauss_copula = seq(0.05, 0.50, by = 0.05),
                 chirp  = seq(0.2, 2.0, by = 0.2),
                 square = seq(0.2, 2.0, by = 0.2),
                 bump   = seq(0.2, 2.0, by = 0.2),
                 trig_copula = seq(0.1, 1.0, by = 0.1))

run_cell <- function(cell, B, CAL) {
  set.seed(SEED + 1009L*B$n + 7919L*cell$ai + 104729L*cell$si)
  g <- gen[[cell$alt]]
  M <- t(vapply(seq_len(R1), function(i) {
    d <- g(B$n, cell$s)
    components(rank(d$x, ties.method = "first"), rank(d$y, ties.method = "first"), B)
  }, numeric(length(comp_names(B$K)))))
  colnames(M) <- comp_names(B$K)
  P <- pvals(M, CAL$laws)
  pw <- vapply(names(CAL$V), function(v) {
    s <- CAL$V[[v]]; ct <- CAL$cal[[v]]$cut
    if (s$type == "single") mean(M[, s$comp] > ct)
    else if (s$type == "sum") mean(rowSums(M[, s$comp]) > ct)
    else mean(do.call(pmin, as.data.frame(P[, s$comp, drop = FALSE])) <= ct)
  }, 0)
  data.frame(n = B$n, alt = cell$alt, strength = cell$s, variant = names(pw),
             power = as.numeric(pw), stringsAsFactors = FALSE)
}

## ------------------------------- run ---------------------------------
## Results go to pilot_results/ beside this script (or beside the working
## directory when run by selection), never into results/, whose cell files
## belong to the main studies.
OUT <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  d <- if (length(m)) dirname(sub("^--file=", "", m[1L])) else getwd()
  file.path(d, "pilot_results")
})
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
res <- list(); nul <- list()
for (n in NS) {
  B <- make_bases(n)
  t0 <- Sys.time(); CAL <- calibrate(B)
  message(sprintf("n = %d: null calibration %.0f s", n, as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  nul[[as.character(n)]] <- data.frame(n = n, variant = names(CAL$cal),
     n_stats = vapply(CAL$V, function(s) length(s$comp), 0L),
     cut = vapply(CAL$cal, `[[`, 0, "cut"), size = vapply(CAL$cal, `[[`, 0, "size"),
     meff = vapply(CAL$cal, `[[`, 0, "meff"), stringsAsFactors = FALSE)
  cells <- list()
  for (ai in seq_along(gen)) for (si in seq_along(STRENGTH[[ai]]))
    cells[[length(cells) + 1L]] <- list(alt = names(gen)[ai], ai = ai, si = si,
                                        s = STRENGTH[[names(gen)[ai]]][si])
  t0 <- Sys.time()
  out <- mclapply(cells, run_cell, B = B, CAL = CAL, mc.cores = NCORE)
  bad <- vapply(out, function(z) !is.data.frame(z), TRUE)
  if (any(bad)) stop("cells failed: ", paste(which(bad), collapse = ","))
  res[[as.character(n)]] <- do.call(rbind, out)
  message(sprintf("n = %d: %d power cells %.0f s", n, length(cells),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  saveRDS(list(CAL = CAL[c("cal", "V", "nullcor")], B = B[c("n", "p", "K")]),
          file.path(OUT, sprintf("pilot_null_n%d.rds", n)))
}
pow <- do.call(rbind, res); nul <- do.call(rbind, nul)
write.csv(pow, file.path(OUT, "pilot_power.csv"), row.names = FALSE)
write.csv(nul, file.path(OUT, "pilot_null.csv"),  row.names = FALSE)
message("done")
