#' Deformation plots: the transformations that maximise the correlation
#'
#' @description
#' `depanalysis()` estimates the pairs of transformations \eqn{(\hat g_s,
#' \hat h_s)} that make \eqn{X} and \eqn{Y} as nearly linearly related as the
#' chosen spans allow, and draws them. The plots of \eqn{(x, \hat g_s(x))}
#' and \eqn{(y, \hat h_s(y))} are the *deformation plots* of Section 2.7 of the
#' paper: they display the reshaping of the two variables that maximises
#' their correlation, and their shape is often enough to recognise the form
#' of the dependence.
#'
#' @details
#' **What is plotted.** Writing \eqn{\tilde g} and \eqn{\tilde h} for the
#' transformations of the copula scale, \eqn{\hat g(\cdot) = \tilde
#' g\{\hat F_X(\cdot)\}} is read off a pair of sample canonical variates
#' with \eqn{F_X} replaced by its empirical counterpart, and likewise for
#' \eqn{\hat h}. The construction is not confined to the leading pair: each
#' subsequent pair is uncorrelated with those before it and carries
#' correlation \eqn{\hat\lambda_s}, so `dimensions = 1:3` draws three modes
#' of association side by side, each panel labelled with its own
#' \eqn{\hat\lambda_s}. The first pair is the natural place to start, since
#' it attains the strongest correlation, but every pair with a
#' non-negligible correlation repays inspection.
#'
#' **Scaling.** The variates are standardised to unit variance, so a panel
#' reads on a scale that does not depend on \eqn{n}, amplitudes are
#' comparable across pairs and across data sets, and the scatter of
#' \eqn{(\hat g_s, \hat h_s)} has slope exactly \eqn{\hat\lambda_s}.
#'
#' **Sign.** A canonical pair is determined only up to a joint change of
#' sign: \eqn{(-\hat g_s, -\hat h_s)} has the same correlation and is the
#' same solution. Left to the singular value decomposition the choice is
#' arbitrary, so the same data can produce a vertically mirrored plot from
#' one run, one order or one implementation to the next. Each pair is
#' therefore oriented so that \eqn{\hat g_s} averages positive over the
#' upper decile of \eqn{x} -- the transform is drawn rising at the
#' right-hand end of the support -- and \eqn{\hat h_s} is flipped alongside
#' it, since flipping one alone would reverse the sign of
#' \eqn{\hat\lambda_s}. No information is lost: the sign of a pair carries
#' none, and what does carry information, the orientation of \eqn{\hat h_s}
#' relative to \eqn{\hat g_s}, is preserved.
#'
#' The rule is a property of the fitted function rather than of its
#' coefficients, so it gives the same answer whichever basis produced the
#' shape. Two obvious alternatives do not: making \eqn{\hat g_s} positively
#' correlated with \eqn{x} is a coin toss whenever \eqn{\hat g_s} is
#' symmetric, which is the V and cosine cases these plots are most used
#' for; and making the observation of largest \eqn{|\hat g_s|} positive
#' turns on which of two extremes happens to be the larger, and was
#' measured flipping between \eqn{p = 5} and \eqn{p = 6} on an otherwise
#' unchanged example. Where the upper decile of \eqn{\hat g_s} is itself
#' near zero -- dependence confined to the left tail -- the rule falls back
#' to making the largest basis coefficient positive. Use `flip` to override
#' either.
#'
#' **Relation to [indeptest()].** The statistic and p-value in the caption
#' come from [indeptest()] itself, called on the same arguments, so they are
#' the numbers that function reports and nothing is recomputed differently
#' here. The canonical correlations in `$cor` are obtained from a Cholesky
#' whitening, which is invariant to the scale and to any non-singular linear
#' reparametrisation of either basis, so with `basis = "poly"` they agree
#' with [indeptest()] to machine precision. With `basis = "spline"` and
#' \eqn{p \neq q} they can differ slightly, because this function builds
#' each margin's basis at its own order while [indeptest()] takes the
#' leading columns of the basis at \eqn{\max(p,q)}; the spline basis is not
#' nested, so those are not the same functions. The per-margin construction
#' is used here because a transformation drawn from the leading few
#' B-splines of a higher-order basis is supported near one end of the range
#' and makes a poor picture.
#'
#' @param x,y Numeric vectors of the same length. Factors are not supported:
#'   a curve through \eqn{(x, \hat g(x))} presupposes an ordering of the
#'   support. Use [indeptest()] for the test itself in that case.
#' @param p,q Orders of the two bases, defaulting to rule (12) of the paper,
#'   \eqn{\max(1, \lfloor n^{0.31}\rfloor - 1)}, with `q = p`.
#' @param basis,ties As in [indeptest()].
#' @param dimensions Integer vector giving the canonical pairs to draw.
#'   Defaults to the first.
#' @param annotate Which test to report in the plot caption: `"Bartlett"`
#'   (the default), `"Pillai"`, `"both"`, or `"none"`. It is placed in the
#'   caption rather than the subtitle so that adding a title with
#'   `ggtitle()` cannot silently displace it.
#' @param flip Logical, recycled over `dimensions`: reverse the sign of a
#'   pair after the convention above has been applied. The sign carries no
#'   information, so this only changes which way up the picture is drawn.
#' @param xlab,ylab Axis labels for the untransformed variables. The
#'   transformed axes are labelled \eqn{g(x)} and \eqn{h(y)}.
#'
#' @return An object of class `"depanalysis"`, a list with
#'   \item{`stat`,`pvalue`}{named vectors holding the Bartlett and Pillai
#'     statistics and their p-values, as returned by [indeptest()].}
#'   \item{`cor`}{all the sample canonical correlations \eqn{\hat\lambda_s},
#'     decreasing. The first is the sample maximal correlation.}
#'   \item{`g`,`h`}{\eqn{n \times} `length(dimensions)` matrices of the
#'     transformed observations, in the order of `x` and `y` as supplied, at
#'     unit variance and with the sign convention above.}
#'   \item{`coef_x`,`coef_y`}{the basis coefficients of the transformations,
#'     for anyone who wants to evaluate them elsewhere.}
#'   \item{`plotx`,`ploty`,`plotxy`}{ggplot objects: the two deformation
#'     plots and the scatter of the transformed pair. `NULL`, with a
#'     warning, if ggplot2 is not installed.}
#'
#' @seealso [indeptest()] for the test, [indepscan()] for the order scan.
#'
#' @examples
#' set.seed(2022)
#' n <- 1000
#'
#' ## X is the conditional standard deviation of Y: a V-shape
#' x <- sqrt(rchisq(n, df = 2)); y <- rnorm(n) * x
#' m <- depanalysis(x, y)
#' m$plotx
#' m$plotxy
#'
#' ## a noisy circle: the transformations are recognisably trigonometric
#' th <- runif(n, 0, 2 * pi)
#' m2 <- depanalysis(cos(th) + rnorm(n, sd = 0.1),
#'                   sin(th) + rnorm(n, sd = 0.1))
#' m2$plotx
#'
#' ## two modes of association at once
#' m2$plotxy
#' depanalysis(cos(th) + rnorm(n, sd = .1), sin(th) + rnorm(n, sd = .1),
#'             p = 6, dimensions = 1:3)$plotx
#' @export
depanalysis <- function(x, y,
                        p = max(1L, as.integer(floor(length(x)^0.31)) - 1L),
                        q = p,
                        basis = c("poly", "spline"),
                        ties  = c("first", "last", "random", "average"),
                        dimensions = 1L,
                        annotate = c("Bartlett", "Pillai", "both", "none"),
                        flip = FALSE,
                        xlab = "x", ylab = "y") {
  basis    <- match.arg(basis)
  ties     <- match.arg(ties)
  annotate <- match.arg(annotate)
  n <- length(x)
  if (length(y) != n)
    stop("x and y must have the same length.", call. = FALSE)
  if (is.factor(x) || is.factor(y) || is.character(x) || is.character(y))
    stop("depanalysis() needs numeric x and y: a curve through (x, g(x)) ",
         "presupposes an ordering of the support. Use indeptest() for the ",
         "test on categorical variables.", call. = FALSE)
  if (!is.numeric(x) || !is.numeric(y))
    stop("x and y must be numeric.", call. = FALSE)
  ok <- stats::complete.cases(x, y)
  if (!all(ok)) {
    warning(sum(!ok), " incomplete case(s) removed.", call. = FALSE)
    x <- x[ok]; y <- y[ok]; n <- length(x)
  }
  force(q)                    # default q = p, taken before p is capped
  p <- cap_order(x, p, "x", !missing(p))
  q <- cap_order(y, q, "y", !missing(q))
  r <- min(p, q)
  dimensions <- as.integer(dimensions)
  if (!length(dimensions) || anyNA(dimensions) ||
      any(dimensions < 1L) || any(dimensions > r))
    stop(sprintf("dimensions must be integers in 1:%d, since min(p, q) = %d.",
                 r, r), call. = FALSE)

  ## ---- the two bases, each at its own order -------------------------
  bf <- if (basis == "poly") chebyshev_basis else ortho_spline_basis_int
  rx <- rank(x, ties.method = ties)
  ry <- rank(y, ties.method = ties)
  U  <- suppressWarnings(bf(n, p))[as.integer(rx), seq_len(p), drop = FALSE]
  V  <- suppressWarnings(bf(n, q))[as.integer(ry), seq_len(q), drop = FALSE]

  ## ---- canonical decomposition, whitened by Cholesky ----------------
  ## Unconditionally, not only when there are ties.  It costs nothing at
  ## these orders and is invariant to the scale of the basis, so this
  ## function cannot be affected by a basis whose normalisation does not
  ## match what the fast path in indeptest() assumes.
  U <- U - rep(colMeans(U), each = n)
  V <- V - rep(colMeans(V), each = n)
  Ru <- tryCatch(chol(crossprod(U) / n), error = function(e)
    stop("the x basis is singular at order p = ", p,
         ": lower p, or use a variable with more distinct values.",
         call. = FALSE))
  Rv <- tryCatch(chol(crossprod(V) / n), error = function(e)
    stop("the y basis is singular at order q = ", q,
         ": lower q, or use a variable with more distinct values.",
         call. = FALSE))
  iRu <- backsolve(Ru, diag(p))
  iRv <- backsolve(Rv, diag(q))
  sv  <- svd(crossprod(iRu, crossprod(U, V) / n) %*% iRv)
  lambda <- pmin(pmax(sv$d, 0), 1)

  A <- iRu %*% sv$u                      # basis coefficients of g_s
  Bc <- iRv %*% sv$v                     # basis coefficients of h_s
  G <- U %*% A                           # unit variance by construction
  H <- V %*% Bc

  ## Sign: see the details section.  g and h must flip together, or the
  ## canonical correlation would change sign.
  hi <- x >= stats::quantile(x, 0.9)
  for (s in seq_len(ncol(G))) {
    tail_mean <- mean(G[hi, s])
    sgn <- if (abs(tail_mean) > 0.1 * stats::sd(G[, s]))
             sign(tail_mean) else sign(A[which.max(abs(A[, s])), s])
    if (sgn < 0) {
      G[, s] <- -G[, s]; H[, s] <- -H[, s]
      A[, s] <- -A[, s]; Bc[, s] <- -Bc[, s]
    }
  }
  if (!identical(flip, FALSE)) {
    fl <- rep_len(as.logical(flip), ncol(G))
    fl[is.na(fl)] <- FALSE
    for (s in which(fl)) {
      G[, s] <- -G[, s]; H[, s] <- -H[, s]
      A[, s] <- -A[, s]; Bc[, s] <- -Bc[, s]
    }
  }
  G <- G[, dimensions, drop = FALSE]; H <- H[, dimensions, drop = FALSE]
  A <- A[, dimensions, drop = FALSE]; Bc <- Bc[, dimensions, drop = FALSE]
  nm <- paste0("dim", dimensions)
  colnames(G) <- colnames(H) <- colnames(A) <- colnames(Bc) <- nm

  ## ---- the test, from indeptest() itself ----------------------------
  tst <- indeptest(x, y, p = p, q = q, basis = basis, ties = ties,
                   test = c("Bartlett", "Pillai"))
  stat <- c(B_n = tst$B_stat, P_n = tst$P_stat)
  pval <- c(B_n = tst$B_pvalue, P_n = tst$P_pvalue)

  out <- list(stat = stat, pvalue = pval, n = n, p = p, q = q,
              basis = basis, ties = ties, cor = lambda,
              dimensions = dimensions, g = G, h = H,
              coef_x = A, coef_y = Bc,
              method = "Pelagatti-Monti deformation analysis")
  pl <- dep_plots(x, y, out, annotate, xlab, ylab)
  structure(c(out, pl), class = "depanalysis")
}


#' The three ggplot objects of a deformation analysis
#'
#' @description
#' Separated from the computation so that ggplot2 stays a suggestion rather
#' than a dependency: the test is the package's reason to exist and must not
#' pull in a plotting stack. Builds its data frames in base R; tidyr is not
#' needed to stack a few columns.
#'
#' @keywords internal
#' @noRd
## The aes() calls below refer to columns of data frames built here, not to
## objects: declared so R CMD check does not read them as undefined globals.
## The .data pronoun would be the modern alternative, but it belongs to
## rlang and reaching it as ggplot2::.data breaks the pronoun inside vars().
utils::globalVariables(c("v", "f", "panel", "gx", "hy"))

dep_plots <- function(x, y, z, annotate, xlab, ylab) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    warning("ggplot2 is not installed, so no plots were made. ",
            "The transformations are in $g and $h.", call. = FALSE)
    return(list(plotx = NULL, ploty = NULL, plotxy = NULL))
  }
  n <- z$n; d <- z$dimensions
  lab <- sprintf("dim%d   cor = %.2f", d, z$cor[d])
  cap <- dep_caption(z, annotate)

  ## geom_line joins points in the order of the data, so x and y are sorted.
  ox <- order(x); oy <- order(y)
  dfx <- data.frame(v   = rep(x[ox], length(d)),
                    panel = factor(rep(lab, each = n), levels = lab),
                    f   = as.vector(z$g[ox, , drop = FALSE]))
  dfy <- data.frame(v   = rep(y[oy], length(d)),
                    panel = factor(rep(lab, each = n), levels = lab),
                    f   = as.vector(z$h[oy, , drop = FALSE]))
  dfxy <- data.frame(gx  = as.vector(z$g),
                     hy  = as.vector(z$h),
                     panel = factor(rep(lab, each = n), levels = lab))

  curve <- function(df, xl, yl)
    ggplot2::ggplot(df, ggplot2::aes(x = v, y = f)) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free_y") +
      ggplot2::labs(x = xl, y = yl, caption = cap)

  list(
    plotx = curve(dfx, xlab, "g(x)"),
    ploty = curve(dfy, ylab, "h(y)"),
    plotxy = ggplot2::ggplot(dfxy, ggplot2::aes(x = gx, y = hy)) +
      ggplot2::geom_point(alpha = 0.4) +
      ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      ggplot2::facet_wrap(ggplot2::vars(panel), scales = "free") +
      ggplot2::labs(x = "g(x)", y = "h(y)", caption = cap))
}


#' Caption text for the deformation plots
#' @keywords internal
#' @noRd
dep_caption <- function(z, annotate) {
  if (identical(annotate, "none")) return(NULL)
  df <- z$p * z$q
  one <- function(k) {
    pv <- format.pval(z$pvalue[[k]], digits = 3, eps = 1e-4)
    sprintf("%s = %.4g on %d df,  %s", k, z$stat[[k]], df,
            if (startsWith(pv, "<")) paste("p", pv) else paste("p =", pv))
  }
  switch(annotate,
         Bartlett = one("B_n"),
         Pillai   = one("P_n"),
         both     = paste(one("B_n"), one("P_n"), sep = ";   "))
}


#' @export
print.depanalysis <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat(sprintf("n = %d,  p = %d,  q = %d,  basis = %s,  ties = %s\n",
              x$n, x$p, x$q, x$basis, x$ties))
  cat(sprintf("\nB_n = %.4g on %d df,  p-value %s\n", x$stat[["B_n"]],
              x$p * x$q, format.pval(x$pvalue[["B_n"]], digits = 3, eps = 1e-4)))
  cat(sprintf("P_n = %.4g on %d df,  p-value %s\n", x$stat[["P_n"]],
              x$p * x$q, format.pval(x$pvalue[["P_n"]], digits = 3, eps = 1e-4)))
  k <- min(length(x$cor), 6L)
  cat("\ncanonical correlations:",
      paste(sprintf("%.3f", x$cor[seq_len(k)]), collapse = "  "),
      if (length(x$cor) > k) " ..." else "", "\n", sep = " ")
  cat(sprintf("\ndrawn: %s  (in $plotx, $ploty, $plotxy)\n",
              paste0("dim", x$dimensions, collapse = ", ")))
  invisible(x)
}
