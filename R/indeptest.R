#' Test of independence between two variables
#'
#' @description
#' Tests the null hypothesis of independence between
#' two random variables \code{x} and \code{y}. The
#' method used depends on the types of the two variables:
#' \itemize{
#'   \item \strong{numeric vs numeric}: orthonormal basis
#'         expansion test based on canonical correlations
#'         between polynomial or B-spline bases applied
#'         to the ranks of \code{x} and \code{y}. Integer
#'         vectors are treated as numeric.
#'   \item \strong{factor vs factor}: same test with
#'         polynomial encoding of factor levels; equivalent
#'         to a chi-square test of association on the
#'         contingency table. Logical and character vectors
#'         are converted to factors.
#'   \item \strong{factor vs numeric} (and the symmetric
#'         case): mixed encoding, with a polynomial basis
#'         for the factor and a spline or polynomial basis
#'         for the numeric variable.
#' }
#'
#' Under the null hypothesis of independence, the Pillai
#' and Bartlett test statistics both converge in
#' distribution to a chi-square with \eqn{pq} degrees of
#' freedom, where \eqn{p} and \eqn{q} are the numbers of
#' basis functions used for \code{x} and \code{y}
#' respectively.
#'
#' @param x A numeric vector or factor.
#' @param y A numeric vector or factor of the same length
#'   as \code{x}.
#' @param p Integer. Number of basis functions for \code{x}.
#'   Defaults depend on the types of the variables. It is
#'   reduced to the number of distinct values of \code{x}
#'   minus one whenever it exceeds that bound, since the
#'   expansion cannot have larger rank; a warning is issued
#'   only if the value was supplied by the user.
#' @param q Integer. Number of basis functions for \code{y}.
#'   Defaults depend on the types of the variables, and are
#'   capped as described for \code{p}.
#' @param basis Character string specifying the type of basis functions
#'   (\code{"poly"}, \code{"spline"}, or \code{"dummy"}, the latter being
#'   valid only for factors).
#' @param basis_fct Character string specifying the type of basis functions
#'   for factor variables (\code{"poly"}, \code{"spline"}, or \code{"dummy"}).
#' @param basis_num Character string specifying the type of basis functions
#'   for numeric variables (\code{"poly"} or \code{"spline"}).
#' @param test Vector of character strings: can be \code{"Pillai"},
#'   \code{"Bartlett"}. Multiple selection is allowed.
#'   Selects the test statistics to compute.
#' @param ties Method to manage ties in numeric variables: choose among
#'   \code{"random"}, \code{"first"}, \code{"last"}, or \code{"average"}.
#'   Corresponds to the \code{ties.method} parameter in the \code{rank()} function.
#'   Only \code{"average"}, or any other convention giving tied observations a
#'   common rank, makes the numeric encoding of a categorical variable
#'   equivalent to its one-hot encoding; the other choices let the arbitrary
#'   ordering within a group of ties enter the statistic.
#' @param ... Additional arguments passed to the appropriate method.
#'
#' @return An object of class \code{"indeptest"}, which
#'   is a list with components:
#'   \describe{
#'     \item{\code{P_stat}}{Value of the Pillai
#'       test statistic \eqn{P_n}.}
#'     \item{\code{B_stat}}{Value of the Bartlett
#'       test statistic \eqn{B_n}.}
#'     \item{\code{P_pvalue}}{P-value for the Pillai
#'       statistic.}
#'     \item{\code{B_pvalue}}{P-value for the Bartlett
#'       statistic.}
#'     \item{\code{method}}{Character string describing
#'       the method used.}
#'     \item{\code{var_types}}{Character vector of length 2
#'       giving the types of \code{x} and \code{y}.}
#'     \item{\code{p}}{Number of basis functions used
#'       for \code{x}.}
#'     \item{\code{q}}{Number of basis functions used
#'       for \code{y}.}
#'     \item{\code{basis}}{String with the name of the basis functions.}
#'     \item{\code{nobs}}{Sample size.}
#'   }
#'
#' @seealso
#'   \code{\link{print.indeptest}} for printing results.
#'
#' @references
#'   Monti, G.S. and Pelagatti, M. (2024). A nonparametric
#'   test of independence between two random variables of
#'   any kind. \emph{Unpublished manuscript}.
#'
#' @examples
#' ## numeric vs numeric
#' set.seed(1)
#' x <- rnorm(200)
#' y <- x^2 + rnorm(200, sd = 0.5)
#' indeptest(x, y)
#'
#' ## factor vs factor
#' x_f <- factor(sample(letters[1:4], 200, replace = TRUE))
#' y_f <- factor(sample(letters[1:3], 200, replace = TRUE))
#' indeptest(x_f, y_f)
#'
#' ## factor vs numeric
#' indeptest(x_f, y)
#'
#' ## integer codes of a nominal variable: with average ranks the numeric
#' ## method reproduces the one-hot encoding exactly
#' x_i <- as.integer(x_f)
#' y_i <- as.integer(y_f)
#' indeptest(x_i, y_i, ties = "average")
#' indeptest(x_f, y_f, basis = "dummy")
#'
#' @importFrom stats pchisq poly
#' @export
indeptest <- function(x, y, ...) {
  if (anyNA(x) || anyNA(y)) stop("Missing values (NA) are not allowed in indeptest.")
  UseMethod("indeptest")
}
#' Cap an expansion order at the number of distinct values minus one
#'
#' @description
#' The scaled ranks of a variable taking \eqn{c} distinct values span a space
#' of dimension at most \eqn{c - 1} once the constant is removed, so an order
#' larger than \eqn{c - 1} makes the cross-product matrix singular when tied
#' observations share a rank, and describes nothing but the arbitrary ordering
#' within groups of ties otherwise.
#'
#' @param v The variable whose distinct values bound the order.
#' @param k The requested order.
#' @param label Name of the variable, used in messages.
#' @param explicit Logical: was the order supplied by the user rather than
#'   taken from the default rule? A warning is issued only if it was.
#'
#' @return The order, reduced if necessary, as a positive integer.
#' @keywords internal
#' @noRd
cap_order <- function(v, k, label = "x", explicit = FALSE) {
  m <- length(unique(v)) - 1L
  if (m < 1L)
    stop(sprintf("%s is constant: the test needs at least two distinct values.",
                 label), call. = FALSE)
  k <- as.integer(k)
  if (k > m) {
    if (explicit)
      warning(sprintf(
        "%s: order reduced from %d to %d, the number of distinct values minus one.",
        label, k, m), call. = FALSE)
    k <- m
  }
  max(1L, k)
}
#' Whitened cross-product of two design matrices
#'
#' @description
#' Returns \eqn{S_{uu}^{-1/2} S_{uv} S_{vv}^{-1/2}}, whose singular values are
#' the sample canonical correlations between the columns of \code{U} and those
#' of \code{V}.
#'
#' When the scaled ranks are a permutation of the grid, the design matrix is a
#' row permutation of the basis and therefore inherits its centring and its
#' orthonormality: \eqn{S_{uu}} and \eqn{S_{vv}} are identity matrices and the
#' cross-product alone is already the required matrix. That is the case for
#' data without ties, and for every tie convention that separates tied
#' observations. It is not the case when tied observations share a rank, as
#' with \code{ties = "average"}: only as many rows of the basis are then used
#' as there are distinct values, each with the multiplicity of its group, so
#' neither property survives and the two moment matrices must be formed and
#' inverted explicitly.
#'
#' @param U,V Design matrices with the same number of rows.
#' @param whiten Logical. If \code{FALSE} the cross-product alone is returned,
#'   which is exact when the ranks are a permutation of the grid and much
#'   cheaper, since it avoids two \eqn{n p^2} products.
#'
#' @return A matrix with \code{ncol(U)} rows and \code{ncol(V)} columns.
#' @keywords internal
#' @noRd
whiten_cross <- function(U, V, whiten) {
  n <- nrow(U)
  if (!whiten)
    return(crossprod(U, V)/n)
  U <- U - rep(colMeans(U), each = n)
  V <- V - rep(colMeans(V), each = n)
  Ru <- chol(crossprod(U)/n)
  Rv <- chol(crossprod(V)/n)
  t(backsolve(Rv, t(backsolve(Ru, crossprod(U, V)/n, transpose = TRUE)),
              transpose = TRUE))
}
#' Internal double-dispatch helper for indeptest
#'
#' @description
#' Builds the name of the second-level method from the implicit classes of the
#' two variables. Integer vectors are mapped to \code{"numeric"}, since the
#' storage mode does not distinguish a count from an ordinal score or a nominal
#' code, and the order cap of \code{cap_order} makes the numeric path correct
#' in all three cases. Logical vectors are converted to factors.
#'
#' @param x,y The two variables passed to \code{indeptest}.
#' @param ... Further arguments.
#' @keywords internal
#' @noRd
indeptest_dispatch2 <- function(x, y, ...) {
  cl_x <- if(is.factor(x)) {
    x <- droplevels(x)
    "factor"
  } else if (is.logical(x)) {
    x <- factor(x)
    "factor"
  } else if (is.character(x)) {
    "character"
  } else if (is.numeric(x)) {
    "numeric"
  } else class(x)[1]
  cl_y <- if(is.factor(y)) {
    y <- droplevels(y)
    "factor"
  } else if (is.logical(y)) {
    y <- factor(y)
    "factor"
  } else if (is.character(y)) {
    "character"
  } else if (is.numeric(y)) {
    "numeric"
  } else class(y)[1]
  method_name <- paste0("indeptest.", cl_x, ".", cl_y)
  method <- tryCatch(
    get(method_name,
        envir = parent.frame(),
        mode  = "function"),
    error = function(e) NULL
  )
  if (!is.null(method))
    return(method(x, y, ...))
  stop(sprintf(
    "No indeptest method for combination (%s, %s).",
    cl_x, cl_y))
}
# -------------------------------------------------------
# First-level S3 methods
# (minimal documentation: just point to the generic page)
# -------------------------------------------------------
# #' @rdname indeptest
#' @export
indeptest.numeric <- function(x, y, ...) {
  indeptest_dispatch2(x, y, ...)
}
# #' @rdname indeptest
#' @export
indeptest.integer <- function(x, y, ...) {
  indeptest_dispatch2(x, y, ...)
}
# #' @rdname indeptest
#' @export
indeptest.logical <- function(x, y, ...) {
  indeptest_dispatch2(x, y, ...)
}
# #' @rdname indeptest
#' @export
indeptest.factor <- function(x, y, ...) {
  indeptest_dispatch2(x, y, ...)
}
# #' @rdname indeptest
#' @export
indeptest.default <- function(x, y, ...) {
  indeptest_dispatch2(x, y, ...)
}
# -------------------------------------------------------
# Second-level methods: one documentation block each,
# all sharing the same @rdname so they appear on one page.
# -------------------------------------------------------
#' @rdname indeptest
#'
#' @section Method: numeric vs numeric:
#' Both \code{x} and \code{y} are numeric, integer vectors being treated as
#' numeric. The test uses an orthonormal polynomial or cubic B-spline basis
#' applied to the ranks of \code{x} and \code{y}. The default basis
#' is polynomial, and the default number of basis
#' functions follows the rule
#' \eqn{p = q = \max(1, \lfloor n^{0.31} \rfloor - 1)}.
#' Each order is then reduced, if necessary, to the number of distinct values
#' of the corresponding variable minus one. A variable with few distinct
#' values is therefore handled correctly without being declared a factor: with
#' \code{ties = "average"} and the resulting order, the polynomial expansion
#' spans the same space as the one-hot encoding and the two give the identical
#' statistic. When tied observations share a rank the design matrices are no
#' longer row permutations of the basis, so they are centred and whitened
#' explicitly before the canonical correlations are computed.
#'
#' @export
indeptest.numeric.numeric <- function(x, y,
                                      p      = max(1L, floor(length(x)^(0.31)) - 1L),
                                      q      = p,
                                      basis  = c("poly", "spline"),
                                      test   = c("Pillai", "Bartlett"),
                                      ties   = c("random", "first", "last", "average"),
                                      ...) {
  basis <- match.arg(basis)
  test  <- match.arg(test, several.ok = TRUE)
  ties  <- match.arg(ties)
  n     <- length(x)
  if (length(y) != n)
    stop("x and y must have the same length.")
  p   <- cap_order(x, p, "x", !missing(p))
  q   <- cap_order(y, q, "y", !missing(q))
  mpq <- max(p, q)
  bf  <- if (basis == "poly") {
    function(n, k) chebyshev_basis(n, k)
  } else {
    function(n, k)
      ortho_spline_basis_int(n, k)
  }
  rx <- rank(x, ties.method = ties)
  ry <- rank(y, ties.method = ties)
  tied <- anyDuplicated(rx) > 0L || anyDuplicated(ry) > 0L
  B  <- bf(n, mpq)
  U  <- B[as.integer(rx), 1:p, drop = FALSE]
  V  <- B[as.integer(ry), 1:q, drop = FALSE]
  S  <- whiten_cross(U, V, tied)
  B_stat <- NA_real_
  B_pval <- NA_real_
  P_stat <- NA_real_
  P_pval <- NA_real_
  if ("Bartlett" %in% test) {
    l2 <- svd(S, 0, 0)$d^2
    B_stat <- (-n + (p+q+3)/2)*sum(log(1-l2))
    B_pval <- pchisq(B_stat, p*q, lower.tail = FALSE)
    if ("Pillai" %in% test) {
      P_stat <- n*sum(l2)
      P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
    }
  } else {
    P_stat <- n*sum(S*S)
    P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
  }
  structure(
    list(P_stat = P_stat,
         B_stat = B_stat,
         P_pvalue = P_pval,
         B_pvalue = B_pval,
         method   = "Pelagatti-Monti independence test",
         var_types = c("numeric", "numeric"),
         p = p,
         q = q,
         basis = basis,
         nobs = n),
    class = "indeptest")
}
#' @rdname indeptest
#'
#' @section Method: factor vs factor:
#' Both \code{x} and \code{y} are factors. Factor levels
#' are encoded as integers if \code{basis = "poly"} or \code{basis = "splines"}
#' and the corresponding basis functions are used. If \code{basis = "dummy"}
#' each factor is converted in a matrix of dummy variables (i.e.
#' one-hot encoding). In the latter case, \eqn{p} and \eqn{q} are set to the
#' number of levels in the respective factors. In the former two cases, \eqn{p}
#' and \eqn{q} must be smaller than the respective number of factor levels. If
#' these constraints are not respected, \eqn{p} and \eqn{q} are set to the
#' number of levels minus one (largest possible value).
#'
#' @export
indeptest.factor.factor <- function(x, y,
                                    p = nlevels(x)-1,
                                    q = nlevels(y)-1,
                                    basis  = c("poly", "spline", "dummy"),
                                    test   = c("Pillai", "Bartlett"), ...) {
  basis <- match.arg(basis)
  test  <- match.arg(test, several.ok = TRUE)
  n <- length(x)
  if (length(y) != n)
    stop("x and y must have the same length.")
  if (nlevels(x) < 2L || nlevels(y) < 2L)
    stop("Both factors must have >= 2 levels.")
  B_stat <- NA_real_
  B_pval <- NA_real_
  P_stat <- NA_real_
  P_pval <- NA_real_
  if (basis == "dummy") {
    Ut <- Matrix::fac2sparse(x)
    Vt <- Matrix::fac2sparse(y)
    p <- nrow(Ut) - 1
    q <- nrow(Vt) - 1
    mu <- Matrix::rowSums(Ut)
    mv <- Matrix::rowSums(Vt)
    UUt <- Ut / sqrt(mu)
    VVt <- Vt / sqrt(mv)
    Suv <- Matrix::tcrossprod(UUt, VVt)
    if ("Bartlett" %in% test) { # Bartlett
      l2 <- svd(Suv, 0, 0)$d[-1]^2
      B_stat <- (-n + (p+q+3)/2)*sum(log(1-l2))
      B_pval <- pchisq(B_stat, p*q, lower.tail = FALSE)
      if ("Pillai" %in% test) { # Bartlett + Pillai
        P_stat <- n*sum(l2)
        P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
      }
    } else { # Pillai only
      P_stat <- n*(sum(Matrix::diag(Matrix::tcrossprod(Suv)))-1)
      P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
    }
  } else {
    xi  <- as.integer(x)
    yi  <- as.integer(y)
    p <- cap_order(xi, p, "x", !missing(p))
    q <- cap_order(yi, q, "y", !missing(q))
    if (basis == "poly") {
      U <- quick_qr(poly(xi, degree = p))
      V <- quick_qr(poly(yi, degree = q))
    } else if (basis == "spline") {
      U <- quick_qr(splines::bs(xi, df = p))
      V <- quick_qr(splines::bs(yi, df = q))
    }
    Suv <- crossprod(U, V)
    if ("Bartlett" %in% test) {
      l2 <- svd(Suv, 0, 0)$d^2
      B_stat <- (-n + (p+q+3)/2)*sum(log(1-l2))
      B_pval <- pchisq(B_stat, p*q, lower.tail = FALSE)
      if ("Pillai" %in% test) { # Bartlett + Pillai
        P_stat <- n*sum(l2)
        P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
      }
    } else { # Pillai only
      P_stat <- n*sum(Suv*Suv)
      P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
    }
  }
  structure(
    list(P_stat = P_stat,
         B_stat = B_stat,
         P_pvalue = P_pval,
         B_pvalue = B_pval,
         method   = "Pelagatti-Monti independence test",
         var_types = c("factor", "factor"),
         p = p,
         q = q,
         basis = basis,
         nobs = n),
    class = "indeptest")
}
#' @rdname indeptest
#'
#' @section Method: character vs character:
#' Both \code{x} and \code{y} are converted into factor and the proper
#' method is called.
#'
#' @export
indeptest.character.character <- function(x, y, ...) {
  indeptest.factor.factor(factor(x), factor(y), ...)
}
#' @rdname indeptest
#'
#' @section Method: factor vs numeric:
#' \code{x} is a factor and \code{y} is numeric. A
#' polynomial basis with \eqn{K_x - 1} functions is used
#' for \code{x}, where \eqn{K_x} is the number of levels
#' of \code{x}. For \code{y}, a polynomial or B-spline
#' basis is used with \code{q} functions as in the
#' numeric vs numeric method. Both orders are reduced, if
#' necessary, to the number of distinct values of the
#' corresponding variable minus one. When \code{y} has ties and these are
#' given a common rank, its design matrix is centred and orthonormalised
#' explicitly, since it is then no longer a row permutation of the basis.
#'
#' @export
indeptest.factor.numeric <- function(x, y,
                                     p = nlevels(x) - 1,
                                     q = max(1L, floor(length(y)^(0.31)) - 1L),
                                     basis_fct = c("poly", "spline", "dummy"),
                                     basis_num = c("poly", "spline"),
                                     test   = c("Pillai", "Bartlett"),
                                     ties   = c("random", "first", "last", "average"),
                                     ...) {
  basis_fct <- match.arg(basis_fct)
  basis_num <- match.arg(basis_num)
  test <- match.arg(test, several.ok = TRUE)
  ties <- match.arg(ties)
  n <- length(x)
  if (length(y) != n)
    stop("x and y must have the same length.")
  if (nlevels(x) < 2L)
    stop("Factor x must have at least 2 levels.")
  q <- cap_order(y, q, "y", !missing(q))
  B_stat <- NA_real_
  B_pval <- NA_real_
  P_stat <- NA_real_
  P_pval <- NA_real_
  # build U for the factor variable
  if (basis_fct == "dummy") {
    Ut <- Matrix::fac2sparse(x)
    su <- sqrt(Matrix::rowSums(Ut))
    U <- Matrix::t(Ut/su)
    p <- nrow(Ut) - 1
  } else {
    xi  <- as.integer(x)
    p   <- cap_order(xi, p, "x", !missing(p))
    if (basis_fct == "poly") {
      U <- quick_qr(poly(xi, degree = p))
    } else if (basis_fct == "spline") {
      U <- quick_qr(splines::bs(xi, df = p))
    }
  }
  # build V for the numeric variable
  bf  <- if (basis_num == "poly") {
    function(n, k) chebyshev_basis(n, k)
  } else {
    function(n, k)
      ortho_spline_basis_int(n, k)
  }
  ry <- rank(y, ties.method = ties)
  B  <- bf(n, q)
  V  <- B[as.integer(ry), 1:q, drop = FALSE]
  if (anyDuplicated(ry) > 0L) {
    V <- V - rep(colMeans(V), each = n)
    V <- t(backsolve(chol(crossprod(V)), t(V), transpose = TRUE))
  } else {
    V <- V/sqrt(n)
  }
  # compute tests
  Suv <- Matrix::as.matrix(Matrix::crossprod(U, V))
  if ("Bartlett" %in% test) {
    l2 <- svd(Suv, 0, 0)$d^2
    B_stat <- (-n + (p+q+3)/2)*sum(log(1-l2))
    B_pval <- pchisq(B_stat, p*q, lower.tail = FALSE)
    if ("Pillai" %in% test) { # Bartlett + Pillai
      P_stat <- n*sum(l2)
      P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
    }
  } else { # Pillai only
    P_stat <- n*sum(Suv*Suv)
    P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
  }
  structure(
    list(P_stat = P_stat,
         B_stat = B_stat,
         P_pvalue = P_pval,
         B_pvalue = B_pval,
         method   = "Pelagatti-Monti independence test",
         var_types = c("factor", "numeric"),
         p = p,
         q = q,
         basis = paste0(basis_fct, "(factor) ", basis_num, "(numeric)"),
         nobs = n),
    class = "indeptest")
}
#' @rdname indeptest
#'
#' @section Method: numeric vs factor:
#' \code{x} is numeric and \code{y} is a factor. This is
#' the symmetric case of the factor vs numeric method.
#' This function inverts the order of x and y, and p and q
#' and calls the method \code{indeptest.factor.numeric()}.
#'
#' @export
indeptest.numeric.factor <- function(x, y,
                                     p = max(1L, floor(length(x)^(0.31)) - 1L),
                                     q = nlevels(y) - 1,
                                     basis_num = c("spline", "poly"),
                                     basis_fct = c("spline", "poly", "dummy"),
                                     test   = c("Pillai", "Bartlett"),
                                     ties   = c("random", "first", "last", "average"),
                                     ...) {
  indeptest.factor.numeric(x = y, y = x, p = q, q = p,
                           basis_fct = basis_fct, basis_num = basis_num,
                           test = test, ties = ties, ...)
}
#' @rdname indeptest
#'
#' @section Method: character vs numeric:
#' \code{x} is converted into a factor and the method for
#' factor and numeric variables is called.
#'
#' @export
indeptest.character.numeric <- function(x, y, ...) {
  indeptest.factor.numeric(factor(x), y, ...)
}
#' @rdname indeptest
#'
#' @section Method: numeric vs. character:
#' \code{y} is converted into a factor and the method for
#' factor and numeric variables is called switching the roles
#' of the two variables.
#'
#' @export
indeptest.numeric.character <- function(x, y, ...) {
  indeptest.factor.numeric(x = factor(y), y = x, ...)
}
# -------------------------------------------------------
# Print method: separate documentation page
# -------------------------------------------------------
#' Print an indeptest object
#'
#' @description
#' Prints a summary of the results of an independence
#' test produced by \code{\link{indeptest}}.
#'
#' @param x   An object of class \code{"indeptest"}.
#' @param ... Currently ignored.
#'
#' @return \code{x} is returned invisibly.
#'
#' @seealso \code{\link{indeptest}}
#'
#' @examples
#' set.seed(1)
#' res <- indeptest(rnorm(100), rnorm(100))
#' print(res)
#'
#' @export
print.indeptest <- function(x, ...) {
  cat("\n", x$method, "\n", sep = "")
  cat(        "  Variable types  :", x$var_types[1], "vs.", x$var_types[2], "\n")
  cat(        "  Basis type      :", x$basis, "\n")
  cat(sprintf("  Basis dimensions: p = %d, q = %d\n",
              x$p, x$q))
  cat(sprintf("  Sample size     : n = %d\n", x$nobs))
  if (!is.na(x$P_stat)) {
    cat(      "  Pillai stat.    :")
    cat(sprintf(" %.4f  (p-value = %.4f)\n",
                x$P_stat, x$P_pvalue))
  }
  if (!is.na(x$B_stat)){
    cat(      "  Bartlett stat.  :")
    cat(sprintf(" %.4f  (p-value = %.4f)\n\n",
                x$B_stat, x$B_pvalue))
  }
  invisible(x)
}
