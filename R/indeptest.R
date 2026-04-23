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
#'         to the ranks of \code{x} and \code{y}.
#'   \item \strong{factor vs factor}: same test with
#'         polynomial encoding of factor levels; equivalent
#'         to a chi-square test of association on the
#'         contingency table.
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
#' @param ... Additional arguments passed to the
#'   appropriate method; see the method-specific
#'   documentation below.
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
#' @export
indeptest <- function(x, y, ...) {
  UseMethod("indeptest")
}


#' Internal double-dispatch helper for indeptest
#'
#' @param x,y The two variables passed to \code{indeptest}.
#' @param ... Further arguments.
#' @keywords internal
#' @noRd
indeptest_dispatch2 <- function(x, y, ...) {
  cl_x        <- class(x)[1]
  cl_y        <- class(y)[1]
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
#' Both \code{x} and \code{y} are numeric. The test uses
#' an orthonormal polynomial or cubic B-spline basis applied to
#' the ranks of \code{x} and \code{y}. The default basis
#' is polynomial, and the default number of basis
#' functions follows the rule
#' \eqn{p = q = \max(1, \lfloor n^{0.3} \rfloor - 1)}.
#'
#' @param p Integer. Number of basis functions for
#'   \code{x}. Defaults to
#'   \code{max(1, floor(n^(0.3)) - 1)}.
#' @param q Integer. Number of basis functions for
#'   \code{y}. Defaults to \code{p}.
#' @param basis Character string, either \code{"poly"}
#'   (default) or \code{"spline"}, specifying the type of
#'   basis functions.
#' @param test Vector of character strings: can be "Pillai",
#'   "Bartlett", or both. It selects the test statistics to compute.
#' @param ties Method to manage ties in x and y: choose among
#'   c("random", "first", "last"), and
#'   compare the parameter ties.method in the function rank()
#'
#' @export
indeptest.numeric.numeric <- function(x, y,
                                      p      = max(1L, floor(length(x)^(0.3)) - 1L),
                                      q      = p,
                                      basis  = c("poly", "spline"),
                                      test   = c("Pillai", "Bartlett"),
                                      ties   = c("random", "first", "last"),
                                      ...) {
  basis <- match.arg(basis)
  test  <- match.arg(test, several.ok = TRUE)
  ties  <- match.arg(ties)
  n     <- length(x)
  mpq <- max(p, q)
  if (length(y) != n)
    stop("x and y must have the same length.")
  bf  <- if (basis == "poly") {
    function(n, k) chebyshev_basis(n, k)
  } else {
    function(n, k)
      ortho_spline_basis_int(n, k)
  }
  rx <- rank(x, ties.method = ties)
  ry <- rank(y, ties.method = ties)
  B  <- bf(n, mpq)
  U  <- B[rx, 1:p]
  V  <- B[ry, 1:q]
  S  <- crossprod(U, V)/n
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
    P_stat <- n*sum(diag(crossprod(S)))
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
#' @param p Integer. Number of basis functions for
#'   \code{x}. Defaults to \code{nlevels(x)-1}.
#' @param q Integer. Number of basis functions for
#'   \code{y}. Defaults to \code{nlevels(y)-1}.
#' @param basis Character string, either \code{"poly"}
#'   (default), \code{"spline"} or \code{"dummy}, specifying the type of
#'   basis functions.
#' @param test Vector of character strings: can be "Pillai",
#'   "Bartlett", or both. It selects the test statistics to compute.
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
  if (p < 1L || q < 1L)
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
    if (p >= nlevels(x)) p <- nlevels(x)-1
    if (q >= nlevels(y)) q <- nlevels(y)-1
    xi  <- as.integer(x)
    yi  <- as.integer(y)
    if (basis == "poly") {
      U <- quick_rq(poly(xi, degree = p))
      V <- quick_rq(poly(yi, degree = q))
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
        P_stat <- n*(sum(diag(tcrossprod(Suv))))
        P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
    }
    if ("Bartlett" %in% test) {
      l2 <- svd(Suv, 0, 0)$d^2
      B_stat <- (-n + (p+q+3)/2)*sum(log(1-l2))
      B_pval <- pchisq(B_stat, p*q, lower.tail = FALSE)
      if ("Pillai" %in% test) {
        P_stat <- n*sum(l2)
        P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
      }
    } else {
      P_stat <- n*sum(diag(crossprod(Suv)))
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
indeptest.character.character <- function(x, y,
                                    p = length(unique(x)) - 1,
                                    q = length(unique(y)) - 1,
                                    basis  = c("poly", "spline", "dummy"),
                                    test   = c("Pillai", "Bartlett"), ...) {
  indeptest.factor.factor(factor(x), factor(y),
                          p = p, q = q, basis = basis, test = test, ...)
}


#' @rdname indeptest
#'
#' @section Method: factor vs numeric:
#' \code{x} is a factor and \code{y} is numeric. A
#' polynomial basis with \eqn{K_x - 1} functions is used
#' for \code{x}, where \eqn{K_x} is the number of levels
#' of \code{x}. For \code{y}, a polynomial or B-spline
#' basis is used with \code{q} functions as in the
#' numeric vs numeric method.
#'
#' @param p Integer. Number of basis functions for
#'   \code{x}. Defaults to \code{nlevels(x)-1}.
#' @param q Integer. Number of basis functions for
#'   \code{y}. Defaults to \code{max(1L, floor(length(y)^(0.3)) - 1L)}.
#' @param basis_fct Character string, either \code{"poly"}
#'   (default), \code{"spline"} or \code{"dummy}, specifying the type of
#'   basis functions for the factor variable.
#' @param basis_num Character string, either \code{"poly"}
#'   (default), \code{"spline"}, specifying the type of
#'   basis functions for the numeric variable.
#' @param test Vector of character strings: can be "Pillai",
#'   "Bartlett", or both. It selects the test statistics to compute.
#' @param ties Method to manage ties in y (numeric): choose among
#'   c("random", "first", "last"), and
#'   compare the parameter ties.method in the function rank()
#'
#' @export
indeptest.factor.numeric <- function(x, y,
                                     p = nlevels(x) - 1,
                                     q = max(1L, floor(length(y)^(0.3)) - 1L),
                                     basis_fct = c("spline", "poly", "dummy"),
                                     basis_num = c("spline", "poly"),
                                     test   = c("Pillai", "Bartlett"),
                                     ties   = c("random", "first", "last"),
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
      ortho_spline_basis_int(k, k)
  }
  ry <- rank(y, ties.method = ties)
  B  <- bf(n, q)
  V  <- B[ry, 1:q]/sqrt(n)
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
    P_stat <- n*(sum(diag(tcrossprod(Suv))))
    P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
  }
  if ("Bartlett" %in% test) {
    l2 <- svd(Suv, 0, 0)$d^2
    B_stat <- (-n + (p+q+3)/2)*sum(log(1-l2))
    B_pval <- pchisq(B_stat, p*q, lower.tail = FALSE)
    if ("Pillai" %in% test) {
      P_stat <- n*sum(l2)
      P_pval <- pchisq(P_stat, p*q, lower.tail = FALSE)
    }
  } else {
    P_stat <- n*sum(diag(crossprod(Suv)))
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
#' @param p Integer. Number of basis functions for
#'   \code{x}. Defaults to \code{max(1L, floor(length(x)^(0.3)) - 1L)}.
#' @param q Integer. Number of basis functions for
#'   \code{y}. Defaults to \code{nlevels(y)-1}.
#' @param basis_num Character string, either \code{"poly"}
#'   (default), \code{"spline"}, specifying the type of
#'   basis functions for the numeric variable.
#' @param basis_fct Character string, either \code{"poly"}
#'   (default), \code{"spline"} or \code{"dummy}, specifying the type of
#'   basis functions for the factor variable.
#' @param test Vector of character strings: can be "Pillai",
#'   "Bartlett", or both. It selects the test statistics to compute.
#' @param ties Method to manage ties in x (numeric): choose among
#'   c("random", "first", "last"), and
#'   compare the parameter ties.method in the function rank()
#'
#' @export
indeptest.numeric.factor <- function(x, y,
                                     p = max(1L, floor(length(x)^(0.3)) - 1L),
                                     q = nlevels(y) - 1,
                                     basis_num = c("spline", "poly"),
                                     basis_fct = c("spline", "poly", "dummy"),
                                     test   = c("Pillai", "Bartlett"),
                                     ties   = c("random", "first", "last"),
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
indeptest.character.numeric <- function(x, y,
                                     p = nlevels(x) - 1,
                                     q = max(1L, floor(length(y)^(0.3)) - 1L),
                                     basis_fct = c("spline", "poly", "dummy"),
                                     basis_num = c("spline", "poly"),
                                     test   = c("Pillai", "Bartlett"),
                                     ties   = c("random", "first", "last"),
                                     ...) {
  indeptest.factor.numeric(factor(x), y,
     p = nlevels(x) - 1,
     q = max(1L, floor(length(y)^(0.3)) - 1L),
     basis_fct = c("spline", "poly", "dummy"),
     basis_num = c("spline", "poly"),
     test   = c("Pillai", "Bartlett"),
     ties   = c("random", "first", "last"),
     ...)
}

#' @rdname indeptest
#'
#' @section Method: numeric vs. character:
#' \code{y} is converted into a factor and the method for
#' factor and numeric variables is called switching the roles
#' of the two variables.
#'
#' @export
indeptest.numeric.character <- function(x, y,
                                        p = nlevels(x) - 1,
                                        q = max(1L, floor(length(y)^(0.3)) - 1L),
                                        basis_fct = c("spline", "poly", "dummy"),
                                        basis_num = c("spline", "poly"),
                                        test   = c("Pillai", "Bartlett"),
                                        ties   = c("random", "first", "last"),
                                        ...) {
  indeptest.factor.numeric(x = factor(y), y = x,
                           p = nlevels(x) - 1,
                           q = max(1L, floor(length(y)^(0.3)) - 1L),
                           basis_fct = c("spline", "poly", "dummy"),
                           basis_num = c("spline", "poly"),
                           test   = c("Pillai", "Bartlett"),
                           ties   = c("random", "first", "last"),
                           ...)
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
  if (!is.null(x$P_stat)) {
    cat(      "  Pillai stat.    :")
    cat(sprintf(" %.4f  (p-value = %.4f)\n",
                x$P_stat, x$P_pval))
  }
  if (!is.null(x$B_stat)){
    cat(      "  Bartlett stat.  :")
    cat(sprintf(" %.4f  (p-value = %.4f)\n\n",
                x$B_stat, x$B_pval))
  }
  invisible(x)
}
