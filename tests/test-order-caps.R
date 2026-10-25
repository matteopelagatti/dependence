## Each order is capped by the distinct values of its own variable.
##
## The default q = p is a promise: evaluated after p had been capped, it made
## the order of y inherit the cap set by x, so that two integer-coded factors
## with 3 and 4 levels were tested at (2, 2) rather than (2, 3), and the
## average-rank statistic no longer reproduced the one-hot one.  Plain
## stopifnot(), so R CMD check runs it without a test dependency.

library(dependence)

set.seed(2)
n  <- 200
f1 <- factor(sample(c("low", "mid", "high"), n, replace = TRUE))
f2 <- factor(sample(c("A", "B", "C", "D"), n, replace = TRUE))

## ------------------------------------------------------------------------
## 1.  Integer codes with average ranks reproduce the one-hot statistic,
##     whichever variable has fewer distinct values.
ref <- indeptest(f1, f2, basis = "dummy")
for (sw in c(FALSE, TRUE)) {
  a <- as.integer(if (sw) f2 else f1)
  b <- as.integer(if (sw) f1 else f2)
  r <- indeptest(a, b, ties = "average")
  stopifnot(r$p == length(unique(a)) - 1L,
            r$q == length(unique(b)) - 1L,
            isTRUE(all.equal(r$P_stat, ref$P_stat)),
            isTRUE(all.equal(r$B_stat, ref$B_stat)))
}

## ------------------------------------------------------------------------
## 2.  depanalysis() follows the same rule.
y <- rnorm(n) + as.integer(f1)
d <- depanalysis(as.integer(f1), y, ties = "average")
stopifnot(d$p == 2L,
          d$q == max(1L, as.integer(floor(n^0.31)) - 1L))
