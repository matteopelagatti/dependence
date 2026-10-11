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
##
## ---------------------------------------------------------------------
##  THE REPRODUCIBILITY CONTRACT, October 2026.  This is the property the
##  whole design now turns on, so it is stated once, precisely.
##
##  Every cell draws its own seed from (n, scheme index, nu, q) and nothing
##  else, and every statistic is evaluated inside with_frozen_rng(), so the
##  generator is at the same state before and after each evaluation block.
##  The sequence of DATASETS a cell sees therefore depends on the cell's
##  seed, on gen[[scheme]], and on the LENGTH of the two loops -- but not on
##  which statistics are being computed.  The loop lengths are R0 (the
##  per-cell null loop, which runs whenever at least one statistic is
##  calibrated per cell) and R1.
##
##  Three things follow, and they are what make the study extensible:
##
##   *  A statistic can be ADDED later and evaluated on exactly the same
##      datasets the others saw, so the comparison is PAIRED rather than
##      merely drawn from the same law.  02b_topup_statistics.R does that,
##      and verifies it by recomputing a statistic already on disk and
##      requiring the stored power back to the last digit.
##
##   *  The same holds for the distribution-free critical values: the null
##      datasets in the critdf loop are generated outside the frozen block
##      and the loop runs max(caps) = R0_DF times whatever DF_NAMES holds,
##      so a critical value computed later comes from the same null sample
##      as the ones already stored.
##
##   *  Order of execution, number of cores and cell scheduling are all
##      irrelevant to the numbers.  The parallel run is reproducible.
##
##  VERIFIED, not assumed: the cells written by the October run reproduce
##  exactly under replay at n = 50 and n = 100 (B_n power 0.8642 and 0.3683,
##  recomputed to the last digit).  Keep it that way -- changing R0, R1,
##  SEED, the seed formula, or any gen[[...]] body breaks the contract and
##  invalidates every cell on disk.
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
USE_SCAN <- TRUE      # the grid version of our test, B_scan
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

## Largest order in the scan, K(n).  It MUST grow with n: the scan is
## consistent only against alternatives with energy somewhere in the spans it
## looks at, so a fixed K would cap consistency at that span and
## Theorem~4 would no longer apply.  The rate should be calibrated as rule
## (12) was; pending that, the grid runs from 1 to two orders above the single
## order rule (12) would have chosen, which is K = 4, 5, 6 at n = 50, 100, 200.
## The empirical evidence so far is only two points -- at n = 50 a grid to 4
## beat a grid to 6 (0.470 against 0.461), at n = 200 the two tied -- which is
## the pattern a growing K predicts but does not pin the rate down.  Run the
## sensitivity block at the end of this script before fixing it in the paper.
K_rule <- function(n) as.integer(floor(n^0.31)) + 1L
SCAN_K_GRID <- c(3L, 4L, 5L, 6L, 8L)   # for the sensitivity block only

## Always returns an INTEGER count.  R0, R0_DF and R1 are written as 5e3,
## 2e4, 1e4, which are doubles, so returning the default unchanged made the
## three vapply(..., 0L, ...) calls below fail: vapply coerces integer up to
## double but never double down to integer.
cap_for <- function(nm, default)
  as.integer(if (is.null(EVAL_CAP[[nm]])) default
             else min(default, as.integer(EVAL_CAP[[nm]])))

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

## ---------------------- three further schemes ------------------------
## Added October 2026, to close gaps the referee round exposed.  Set
## EXTRA_SCHEMES to FALSE to reproduce the original ten exactly.  Each adds
## 60 cells to the 600, so all three cost about thirty per cent more.
##
##  lattice        The clustered alternative in the BEAST paper's own
##                 bivariate design is a checkerboard, and it is the one case
##                 there in which BEAST is NOT the top performer, the
##                 chi-squared test being strongest.  Nothing among the ten
##                 has that structure.  The version here is a balanced 4 x 4
##                 lattice, two occupied cells in every row and column,
##                 arranged so that E[j | i] is the same for every i: margins
##                 exactly uniform, and Pearson, Spearman and Kendall
##                 correlation all zero, verified numerically at 4 x 10^5
##                 draws (0.002, 0.002, -0.009).  Every rank-correlation
##                 method is therefore powerless against it by construction.
##
##                 NOTE, because it affects how the result should be read:
##                 the naive checkerboard, mass where i + j is even, is NOT
##                 uncorrelated -- it has Spearman 0.19 at this resolution --
##                 which is why the occupied cells are listed explicitly.
##                 And this lattice is NOT a case where a binary basis should
##                 beat a smooth one: its dependence runs in the
##                 inner-versus-outer direction, and corr{(x-1/2)^2,
##                 (y-1/2)^2} = -0.70, so the quadratic basis function
##                 captures it almost completely and even p = q = 2 should do
##                 well.  Expect to win here.  For a lattice that genuinely
##                 defeats a low-order smooth basis the resolution has to
##                 rise, which needs care to keep the rank correlations at
##                 zero; the high-frequency sine below covers that regime in
##                 the functional direction instead.
##
##  sine 8pi       The existing sine scheme is sin(4 pi x), two cycles.  A
##                 low-order smooth basis is blind to high-frequency
##                 departures by construction -- Proposition 1 says the
##                 noncentrality is the captured energy, which is near zero
##                 when the direction lies outside the span -- and this is
##                 precisely the regime BET, BEAST and Chatterjee's xi are
##                 built for.  Expect to lose here.  That is the point:
##                 after a rejection for not establishing when the test
##                 should be preferred, a design that omits its known worst
##                 case reads as evasion, whereas losing and explaining why,
##                 with the remedy (raise p), is the theory doing its job.
##
##  Clayton        Monotone dependence concentrated in the lower tail.  None
##                 of the ten has tail dependence: "linear" is the only
##                 monotone scheme and its copula is light-tailed.  Sampled
##                 by conditional inversion at Kendall's tau = 1/3, then
##                 given the same additive noise as the rest so that the
##                 noise index keeps its meaning across schemes.
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
## Pseudo-observations, in the same convention as indeptest: ties broken by
## order of appearance, divisor n + 1 so the values stay inside (0, 1).
pseudo <- function(x) rank(x, ties.method = "first") / (length(x) + 1)

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

  ## ---- the grid version of our own test -----------------------------
  ## P_n(k,k) for k = 1..K, each turned into its EXACT null p-value using
  ## the laws tabulated in SCAN_LAW, then the smallest p-value taken.  The
  ## reported quantity is its negative so that, like every other adapter
  ## here, larger means more dependence.
  ##
  ## Equivalently, and this is the form for the paper: reject when
  ## p_k < c w_k for some k, with uniform weights w_k = 1/K and c fixed by
  ## the usual calibration so the size is exactly alpha.  The K statistics
  ## are nested partial sums of one cross-moment array, hence strongly
  ## positively correlated -- null correlations ran from 0.41 to 0.81 at
  ## n = 100 -- so the exact joint calibration is about 1.6 times less
  ## conservative than Bonferroni, and that margin is free power.
  ##
  ## ONE PASS, not K calls to indeptest.  With the polynomial basis and
  ## separated ties the whole profile is a sequence of partial sums of the
  ## squared entries of a single K x K cross-moment array, which is what
  ## dependence:::scan_profile() computes; the basis matrices depend on
  ## (n, K) alone and are built once per n in SCAN_BAS.  Measured cost: 0.41,
  ## 0.50 and 0.57 ms at n = 50, 100, 200 -- the same as ONE evaluation of
  ## B_n at the order rule (12) chooses, because the Bartlett SVD costs about
  ## what the whole profile does.  Looping indeptest instead was 6.3 times
  ## slower and, for basis = "spline", would also have been wrong: that
  ## basis is not nested, since bs() moves its knots with df.  See
  ## ?indepscan and tests/test-indepscan.R in the package.
  ##
  ## scan_profile() is internal, hence the ::: .  It is used rather than the
  ## exported indepscan() because the harness wants the bare statistic and
  ## supplies its own calibration below; indepscan() additionally builds the
  ## weight sweep and the result object, which would be paid 7.8 million
  ## times here for nothing.
  if (USE_SCAN) {
    s$B_scan <- function(x, y, p, q) {
      n <- length(x); key <- as.character(n)
      law <- SCAN_LAW[[key]]
      if (is.null(law)) return(NA_real_)       # before SCAN_LAW is built
      K <- length(law)
      v <- dependence:::scan_profile(x, y, K, "poly", "first", SCAN_BAS[[key]])
      pk <- vapply(seq_len(K), function(k)
        1 - findInterval(v[k], law[[k]]) / length(law[[k]]), 0)
      -min(pk)
    }
  }

  if (USE_NEW) {
    ## BET, binary expansion testing (Zhang, JASA 2019).  The closest
    ## competitor in design: like ours it projects the copula onto a finite
    ## basis and tests the projection, but the basis is binary rather than
    ## smooth, and the aggregation is a maximum over interactions rather
    ## than a sum over a chosen span.  The comparison reads as binary versus
    ## smooth basis at matched resolution, which is the honest framing.
    ##
    ## The package exports no function called BET: the fixed-depth test is
    ## MaxBET, whose depth argument is dep, and the statistic is the extreme
    ## asymmetry.  Pseudo-observations are passed explicitly with
    ## unif.margin = TRUE rather than letting the package apply its own
    ## empirical distribution transform, for three reasons: it is the same
    ## rank convention as indeptest, ties.method = "first" and divisor n + 1,
    ## so the two tests see exactly the same input; it makes the statistic
    ## rank-invariant by construction rather than by assumption; and it
    ## sidesteps the margin handling that BEAST stumbles over below.
    ## THE ABSOLUTE VALUE IS THE FIX, and it is not cosmetic.  Without it the
    ## October run gave BET a rejection rate of 0.003 on the circle, 0.001 on
    ## the lattice and 0.014 on sine 8pi, against a nominal 0.05, with 37 per
    ## cent of its cells below 0.03; on the circle its power ROSE from 0.003
    ## to 0.012 as the dependence weakened.  That is the signature of a signed
    ## quantity calibrated on one tail: on the symmetric schemes the extreme
    ## interaction is negative, so the upper-tail critical value is never
    ## approached.  Extreme.Asymmetry is a signed maximum; the test is on its
    ## magnitude.  If abs() is not enough, z.statistic is the next thing to
    ## try, and MaxBETs() after that -- the orientation check below now
    ## exercises a symmetric alternative, so a wrong choice is caught at
    ## startup rather than after a day of computing.
    s$BET <- function(x, y, p, q)
      abs(pluck1(BET::MaxBET(cbind(pseudo(x), pseudo(y)), dep = d_rule(length(x)),
                             unif.margin = TRUE, asymptotic = TRUE),
                 c("Extreme.Asymmetry", "z.statistic")))

    ## BEAST, the adaptive-weight version, designed for near-uniform power
    ## and the most likely of the four to beat us somewhere.  If it does,
    ## that is the result: say so, and let the calibration and cost
    ## comparison carry the argument instead.  It resamples internally to
    ## build its weights, so its value is random given the data; that extra
    ## variability belongs to the method, not to this harness, and it is why
    ## detect_rank_invariant() sets a common seed before each of its calls.
    ## method = "stat" returns the statistic without building the null
    ## distribution, which is what this harness calibrates for itself; it
    ## skips the B = 100 permutations and is the difference between BEAST
    ## being affordable and not.  The subsampling that forms the adaptive
    ## weights remains, so the value is still random given the data.
    s$BEAST <- function(x, y, p, q)
      pluck1(BET::BEAST(cbind(pseudo(x), pseudo(y)), dep = d_rule(length(x)),
                        unif.margin = TRUE, method = "stat"),
             c("BEAST.Statistic", "Statistic", "statistic"))

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
##
## THE DEFAULT n MATTERS, and it must lie in N_GRID.  B_scan returns NA
## unless SCAN_LAW holds an entry for the n it is handed, and an NA fails
## the comparison below, so the scan would be declared value-based.  That is
## not a cosmetic misclassification: it would be calibrated per cell against
## permuted-y nulls while its per-order laws come from uniforms, and -- far
## worse -- it would join CELL_NAMES, which changes max(caps) in the
## per-cell null loop, which changes the DATASETS every statistic sees.  The
## whole reproducibility contract at the top of this file would break, and
## nothing would say so.  Hard-coding 100 happened to work only because 100
## is in N_GRID; it is now taken from N_GRID, and asserted below.
detect_rank_invariant <- function(n = max(N_GRID), p = p_rule(n), q = p,
                                  reps = 3L, tol = 1e-8) {
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

## Orientation check.  Every adapter must return a LARGER value under strong
## dependence than under independence.  An interface that yields a p-value
## instead of a statistic, or a statistic oriented the other way, passes
## self_test() -- it returns one finite number -- and then silently produces
## near-zero power for the rest of the study, which looks like a result
## rather than a bug.  This is the check that catches it.  pluck1() makes the
## risk real, since it will happily return a p-value field if the statistic
## field is absent under the installed version.
## TWO alternatives, and a statistic must be larger under BOTH.  The October
## run used only the monotone one, which is exactly the case BET handles
## correctly, so a statistic that is reversed on symmetric dependence sailed
## through and then produced power below its own nominal size on a third of
## the design.  The circle is the cheapest symmetric alternative to hand.
orientation_check <- function(n = 200L, p = 3L, q = 3L, reps = 20L) {
  S <- make_stats(); nm <- names(S)
  alts <- list(
    monotone = function(n) { x <- runif(n); list(x = x, y = x + 0.1*rnorm(n)) },
    symmetric = function(n) { x <- runif(n); r <- rbinom(n, 1, 0.5)
      list(x = x, y = (2*r-1)*sqrt(pmax(0, 1-(2*x-1)^2)) + 0.1*rnorm(n)) })
  md <- mi <- matrix(NA_real_, length(alts), length(S),
                     dimnames = list(names(alts), nm))
  for (a in names(alts)) {
    dep <- ind <- matrix(NA_real_, reps, length(S), dimnames = list(NULL, nm))
    for (r in seq_len(reps)) {
      set.seed(5000L + r)
      d <- alts[[a]](n)
      xi <- runif(n); yi <- runif(n)
      for (k in nm) {
        dep[r, k] <- tryCatch(as.numeric(S[[k]](d$x, d$y, p, q)), error = function(e) NA_real_)
        ind[r, k] <- tryCatch(as.numeric(S[[k]](xi, yi, p, q)), error = function(e) NA_real_)
      }
    }
    md[a, ] <- colMeans(dep, na.rm = TRUE); mi[a, ] <- colMeans(ind, na.rm = TRUE)
  }
  ok <- apply(is.finite(md) & is.finite(mi) & md > mi, 2L, all)
  for (k in nm) {
    bad <- names(alts)[!(is.finite(md[, k]) & is.finite(mi[, k]) & md[, k] > mi[, k])]
    message(sprintf("  %-10s monotone %-11.5g / %-11.5g  symmetric %-11.5g / %-11.5g  %s",
      k, md["monotone", k], mi["monotone", k], md["symmetric", k], mi["symmetric", k],
      if (isTRUE(ok[k])) "ok"
      else paste("WRONG ORIENTATION under:", paste(bad, collapse = ", "),
                 "-- a signed statistic, or a p-value, not a magnitude")))
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

## ------------- null laws for the scan, one pass per n ----------------
## THREE independent null samples are involved in the scan and they must
## stay apart:
##
##   (i)   the per-order laws of P_n(k,k), which turn each order's statistic
##         into a p-value -- indepscan_null()'s first batch, $cells;
##   (ii)  indepscan_null()'s second batch, $pmat, which the package uses to
##         calibrate the scan.  The harness does NOT use it for that, since
##         B_scan must be size-adjusted in the same critdf machinery as every
##         other rank-invariant statistic here, on the same null datasets.
##         It is used below only as a free cross-check of (iii);
##   (iii) the critdf null sample at seed SEED + 31 n + 7 p + q, which fixes
##         B_scan's own 5 per cent point alongside B_n's.
##
## Using one sample for two of these would evaluate the scan on the very
## draws that defined its reference distribution, and it would over-reject.
## ------------------ a cell's identity: seed and file name -------------
## Both live ABOVE the sentinel because 02b_topup_statistics.R needs them
## and must not keep its own copy.  A duplicated seed formula is the one
## way the two scripts can disagree without any check noticing: the top-up
## would faithfully reproduce the data of the OLD formula, report success,
## and leave the main script computing something else.
##
## nu enters through round(nu*10) rather than nu itself, so that the seed is
## an exact integer function of the noise index and does not depend on how
## seq() rounded the grid.
cell_seed <- function(n, si, nu, q)
  SEED + 1009L*as.integer(n) + 7919L*as.integer(si) +
    104729L*as.integer(round(nu*10)) + 13L*as.integer(q)

cell_file <- function(cell) file.path(OUTDIR,
  sprintf("cell_n%d_s%02d_nu%03d_q%d.rds", cell$n, cell$si, round(cell$nu*100), cell$q))

## ------------------ fingerprint of a cell's data stream ---------------
## A rejection rate out of R1 is a coarse summary: a change to a generator
## small enough to flip no rejection at all would pass a check based on it.
## These five accumulated moments pin the datasets down directly, and are
## carried as an ATTRIBUTE of the cell -- rbind() drops attributes, so the
## assembled CSV is byte-identical to before -- which lets a later top-up
## verify its replay against the data rather than against a count.  Cells
## written before this was added carry no fingerprint and fall back to the
## rejection-count check, which is what the October cells do.
stream_acc <- function(fp, x, y)
  fp + c(sum(x), sum(y), sum(x*y), sum(x*x), sum(y*y))
stream_match <- function(a, b, tol = 1e-10)
  !is.null(a) && !is.null(b) && length(a) == 5L && length(b) == 5L &&
    all(abs(a - b) <= tol * pmax(1, abs(a), abs(b)))

## A function, not inline code, because 02b_topup_statistics.R needs exactly
## these tables and must not have its own copy of the construction.
build_scan_tables <- function(ns, R = R0_DF) {
  law <- bas <- pmat <- list()
  for (n in ns) {
    K <- K_rule(n); key <- as.character(n)
    ## Basis matrices for orders 1..K: a function of (n, K) alone, so one
    ## build serves the whole study and is handed to every scan_profile()
    ## call rather than being recomputed 7.8 million times.
    bas[[key]] <- dependence:::scan_basis(n, K, "poly")
    f <- file.path(OUTDIR, sprintf("scanlaw_n%d_K%d.rds", n, K))
    if (file.exists(f)) {
      tb <- readRDS(f)
      message(sprintf("scan null laws for n = %d, K = %d: cached", n, K))
    } else {
      message(sprintf("tabulating scan null laws, n = %d, K = %d, R = %g ...",
                      n, K, R))
      t0 <- Sys.time()
      tb <- indepscan_null(n = n, K = K, R = R, basis = "poly",
                           ties = "first", seed = SEED + 101L * n)
      saveRDS(tb, f)
      message(sprintf("  done in %.0f s",
                      as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    }
    stopifnot("cached scan table has the wrong K" = length(tb$cells) == K)
    law[[key]] <- tb$cells; pmat[[key]] <- tb$pmat
  }
  list(law = law, bas = bas, pmat = pmat)
}

## =====================================================================
##  <<< DEFINITIONS END >>>
##
##  Everything ABOVE this line is configuration and function definitions,
##  with no side effect beyond dir.create(OUTDIR).  Everything BELOW it
##  computes.  02b_topup_statistics.R evaluates exactly this prefix, so the
##  two scripts cannot drift apart in the generators, the seed formula or
##  the adapters -- which is what makes a later statistic comparable with
##  the ones already on disk.  Do not move the line below anything that
##  computes, and do not define a generator or an adapter below it.
## =====================================================================

SCAN_LAW <- list(); SCAN_BAS <- list(); SCAN_PMAT <- list()
if (USE_SCAN) {
  .st <- build_scan_tables(N_GRID)
  SCAN_LAW <- .st$law; SCAN_BAS <- .st$bas; SCAN_PMAT <- .st$pmat; rm(.st)
}

## ------------------------- startup diagnostics -----------------------
message("adapter self-test at n = ", N_GRID[1], ":")
STAT_OK <- self_test(N_GRID[1], p_rule(N_GRID[1]), p_rule(N_GRID[1]))
SECS    <- attr(STAT_OK, "secs")
if (!all(STAT_OK)) warning("dropping: ", paste(names(STAT_OK)[!STAT_OK], collapse = ", "))
KEEP <- names(STAT_OK)[STAT_OK]
SECS <- SECS[KEEP]
stopifnot("B_n must be available" = "B_n" %in% KEEP)

message("\norientation check (larger must mean more dependence):")
## at the LARGEST n in the grid, so that B_scan finds its null laws and so
## that the symmetric alternative has enough data to be detectable at all
ORIENT <- orientation_check(n = max(N_GRID), p = p_rule(max(N_GRID)),
                            q = p_rule(max(N_GRID)))[KEEP]
if (!all(ORIENT)) {
  warning("dropping wrongly oriented adapters: ",
          paste(KEEP[!ORIENT], collapse = ", "),
          ". Either the accessor is returning a p-value, in which case negate it, ",
          "or the candidate list in pluck1() is picking the wrong field.")
  KEEP <- KEEP[ORIENT]; SECS <- SECS[KEEP]
}
stopifnot("B_n must be available" = "B_n" %in% KEEP)

message("\nrank-invariance detection (distribution-free under H0):")
RI <- detect_rank_invariant()[KEEP]
for (nm in KEEP) message(sprintf("  %-10s %s", nm,
  if (RI[nm]) "rank-invariant: calibrated once per (n, p, q)" else "value-based: calibrated per cell"))
DF_NAMES  <- names(RI)[RI]; CELL_NAMES <- names(RI)[!RI]

## The classification is load-bearing twice over: it decides how a statistic
## is calibrated, and -- through length(CELL_NAMES) and max(caps) -- how many
## draws the per-cell null loop makes, hence which datasets the alternative
## loop sees.  Two assertions, so that a change here is deliberate.
if ("B_scan" %in% KEEP) stopifnot(
  "B_scan must be rank-invariant; see the note on detect_rank_invariant()" =
    "B_scan" %in% DF_NAMES)
message(sprintf("  per-cell calibration applies to %d statistic(s): %s",
                length(CELL_NAMES),
                if (length(CELL_NAMES)) paste(CELL_NAMES, collapse = ", ") else "none"))
if (length(list.files(OUTDIR, "^cell_.*rds$"))) {
  ex <- tryCatch(readRDS(file.path(OUTDIR, list.files(OUTDIR, "^cell_.*rds$")[1])),
                 error = function(e) NULL)
  if (!is.null(ex)) {
    ## The replication counts stored in the cells are the ones the cells were
    ## actually built with.  If the configuration above has moved since, the
    ## cells will be reported as done and the published table will carry the
    ## old numbers under the new R0/R1 in its caption.
    if (max(ex$R1) != R1)
      warning(sprintf(paste("the cells on disk were built with R1 = %g but R1 is now %g.",
                            "They will be reported as done and NOT recomputed, so the",
                            "table would describe replications the numbers do not have.",
                            "Either restore R1 = %g or delete %s/cell_*.rds."),
                      max(ex$R1), R1, max(ex$R1), OUTDIR), call. = FALSE)
    if (any(!ex$df_calibrated) && max(ex$R0[!ex$df_calibrated]) != R0)
      warning(sprintf(paste("the cells on disk were built with R0 = %g but R0 is now %g.",
                            "R0 sets the length of the per-cell null loop, so it also",
                            "fixes which datasets the alternative loop sees: the cells",
                            "cannot be extended under the new value. Restore R0 = %g or",
                            "delete %s/cell_*.rds."),
                      max(ex$R0[!ex$df_calibrated]), R0,
                      max(ex$R0[!ex$df_calibrated]), OUTDIR), call. = FALSE)
    old <- sort(ex$test[!ex$df_calibrated])
    if (!identical(old, sort(CELL_NAMES)))
      warning(sprintf(paste("the cells on disk were written with %s calibrated per cell,",
                            "now it is %s. The per-cell null loop length changes with that set,",
                            "so the alternative datasets will NOT match the stored cells.",
                            "Delete %s/cell_*.rds and re-run from scratch, or restore the",
                            "previous set of statistics."),
                      paste(old, collapse = ", "), paste(sort(CELL_NAMES), collapse = ", "),
                      OUTDIR), call. = FALSE)
  }
}

## ------------------------- cost projection ---------------------------
PQ <- unique(do.call(rbind, lapply(N_GRID, function(n)
  do.call(rbind, lapply(unique(c(p_rule(n), if (SECONDARY_Q1) 1L)),
                        function(q) c(n, p_rule(n), q))))))
## Timed at EVERY sample size, not only the smallest.  BET and BEAST grow
## with n AND with the binary depth, and d_rule raises the depth from 2 to 3
## between n = 50 and n = 100, multiplying the number of interactions by four
## on top of the growth in n.  A projection from the n = 50 timing alone
## understates them by an order of magnitude, which is exactly the regime in
## which one discovers the problem three days into the run.
time_adapters <- function(n, reps = 5L) {
  S <- make_stats()[KEEP]; p <- p_rule(n)
  set.seed(4242L)
  vapply(KEEP, function(k) {
    t0 <- Sys.time()
    for (r in seq_len(reps)) {
      x <- runif(n); y <- x + 0.5*rnorm(n)
      try(S[[k]](x, y, p, p), silent = TRUE)
    }
    as.numeric(difftime(Sys.time(), t0, units = "secs")) / reps
  }, 0)
}
message("\nseconds per evaluation, at each sample size:")
TIM <- sapply(N_GRID, time_adapters)
if (is.null(dim(TIM))) TIM <- matrix(TIM, nrow = length(KEEP), dimnames = list(KEEP, NULL))
colnames(TIM) <- paste0("n=", N_GRID)
print(signif(TIM, 3))

CELLS_PER_N <- length(gen) * length(NU) *
  (1 + as.integer(SECONDARY_Q1) * length(FUNCTIONAL) / length(gen))
proj <- data.frame(test = KEEP, stringsAsFactors = FALSE)
proj$hours <- vapply(KEEP, function(t)
  sum(vapply(seq_along(N_GRID), function(j) {
    nq <- length(unique(c(p_rule(N_GRID[j]), if (SECONDARY_Q1) 1L)))
    nulls <- if (RI[t]) nq * cap_for(t, R0_DF) else CELLS_PER_N * cap_for(t, R0)
    TIM[t, j] * (CELLS_PER_N * cap_for(t, R1) + nulls)
  }, 0)) / 3600, 0)
proj$share <- proj$hours / sum(proj$hours)
proj <- proj[order(-proj$hours), ]
message("\nprojected cost, from the per-n timings above:")
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

## Cross-check of the scan's calibration, free because both pieces exist.
## The critdf table above fixed B_scan's 5 per cent point from one null
## sample; indepscan_null()'s second batch is an independent one.  Two
## independent estimates of the same quantile should agree to Monte Carlo
## error.  A gap much larger than that means the per-order laws and the
## statistic being calibrated have drifted apart -- a different K, a
## different basis, a different tie convention -- which is the one mistake
## in this construction that would not show up as anything else.
if (USE_SCAN && "B_scan" %in% DF_NAMES) {
  for (n in N_GRID) {
    key <- as.character(n); p <- p_rule(n)
    cv <- crit_df_table[[sprintf("n%d_p%d_q%d", n, p, p)]][["B_scan"]]
    alt <- quantile(-apply(SCAN_PMAT[[key]], 1L, min), 1 - ALPHA, names = FALSE)
    message(sprintf("  B_scan 5%% point, n = %3d: critdf %.5f   independent batch %.5f   gap %.5f",
                    n, cv, alt, abs(cv - alt)))
  }
}

## ------------------------------ one cell -----------------------------
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
  set.seed(cell_seed(n, cell$si, nu, q))

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
  fp <- numeric(5)
  for (i in seq_len(max(caps1))) {
    d <- g(n, nu)
    fp <- stream_acc(fp, d$x, d$y)
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
  attr(out, "stream") <- fp
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
  clusterExport(cl, c("gen","make_stats","run_cell","cell_file","cell_seed",
                      "cell_done","pluck1",
                      "pseudo","with_frozen_rng","N_GRID","R0","R0_DF","R1","ALPHA","SEED",
                      "USE_MIC","USE_NEW","USE_SCAN","p_rule","d_rule","K_rule",
                      "SCAN_LAW","SCAN_BAS","cap_for","EVAL_CAP",
                      "stream_acc","stream_match",
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
  SCALE_TESTS <- setdiff(KEEP, "B_scan")       # see note below
  S <- make_stats()[SCALE_TESTS]
  ns <- c(1e3, 1e4, 1e5, 1e6)
  sc <- expand.grid(test = SCALE_TESTS, n = ns, secs = NA_real_,
                    stringsAsFactors = FALSE)
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
## B_scan is left out of the scaling table deliberately.  Its cost is K(n)
## times that of a single B_n, and its null laws exist only for the n of the
## study, so the figures would be neither meaningful nor comparable.  Quote
## B_n for cost and say in the text that the scan multiplies it by K(n),
## which is 4 to 6 here -- still far below any competitor needing resampling.

## =====================================================================
##  Calibration of K(n), the top of the scan
##
##  K must grow with n: the scan inherits consistency from Theorem 4 only
##  through the spans it actually looks at, so a fixed K caps the set of
##  alternatives it is consistent against.  What the rate should be is an
##  empirical question of the same kind rule (12) answered for a single
##  order, and this block is the study that answers it.  Set it to TRUE and
##  run once; it is not part of the main grid.
##
##  All K variants share one pass: the profile is computed to Kmax and every
##  smaller K is a prefix of the same p-value vector, so the cost is that of
##  the largest K alone, not the sum over the grid.
## =====================================================================
RUN_SCAN_SENS <- FALSE
if (RUN_SCAN_SENS) {
  KMAX <- max(SCAN_K_GRID); R_S <- 2000L; R0_S <- 10000L
  NU_S <- seq(0.2, 2.0, by = 0.2)
  ## One pass to KMAX, with the basis built once; every smaller K is a
  ## prefix of the same vector, so the whole grid costs what the largest K
  ## alone costs.  Same function the main adapter uses.
  BAS_S <- lapply(N_GRID, function(n) dependence:::scan_basis(n, KMAX, "poly"))
  names(BAS_S) <- as.character(N_GRID)
  prof <- function(x, y, K, n = length(x))
    dependence:::scan_profile(x, y, K, "poly", "first", BAS_S[[as.character(n)]])
  sens <- list()
  for (n in N_GRID) {
    f <- file.path(OUTDIR, sprintf("scanlaw_sens_n%d_K%d.rds", n, KMAX))
    if (file.exists(f)) { law <- readRDS(f) } else {
      message(sprintf("sensitivity: null laws n = %d to K = %d", n, KMAX))
      set.seed(SEED + 211L*n)
      M <- t(vapply(seq_len(R0_S), function(i) prof(runif(n), runif(n), KMAX),
                    numeric(KMAX)))
      law <- lapply(seq_len(KMAX), function(k) sort(M[, k]))
      saveRDS(law, f)
    }
    pvec <- function(v) vapply(seq_along(v), function(k)
      1 - findInterval(v[k], law[[k]]) / length(law[[k]]), 0)
    ## critical value for each K, on an independent null batch
    set.seed(SEED + 307L*n)
    P0 <- t(vapply(seq_len(R0_S), function(i) pvec(prof(runif(n), runif(n), KMAX)),
                   numeric(KMAX)))
    cK <- vapply(SCAN_K_GRID, function(K)
      quantile(-apply(P0[, seq_len(K), drop = FALSE], 1L, min), 1 - ALPHA,
               names = FALSE), 0)
    acc <- matrix(0, length(SCAN_K_GRID), 0)
    for (s in names(gen)) for (nu in NU_S) {
      set.seed(SEED + 409L*n + 11L*which(names(gen) == s) + as.integer(nu*10))
      P1 <- t(vapply(seq_len(R_S), function(i) {
        d <- gen[[s]](n, nu); pvec(prof(d$x, d$y, KMAX)) }, numeric(KMAX)))
      acc <- cbind(acc, vapply(seq_along(SCAN_K_GRID), function(j)
        mean(-apply(P1[, seq_len(SCAN_K_GRID[j]), drop = FALSE], 1L, min) > cK[j]), 0))
    }
    sens[[as.character(n)]] <- rowMeans(acc)
    message(sprintf("n = %3d   " , n),
            paste(sprintf("K=%d: %.3f", SCAN_K_GRID, rowMeans(acc)), collapse = "   "))
  }
  out <- do.call(rbind, lapply(names(sens), function(n)
    data.frame(n = as.integer(n), K = SCAN_K_GRID, mean_power = sens[[n]])))
  write.csv(out, paste0(OUTSTEM, "_scanK.csv"), row.names = FALSE)
  message("Written: ", OUTSTEM, "_scanK.csv -- pick K(n) from this, ",
          "then set K_rule accordingly and rerun the main grid.")
}
