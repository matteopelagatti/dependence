src <- readLines("basis_comparison.R"); cut <- grep("^res <- list\\(\\); siz <- list\\(\\)", src)
eval(parse(text = src[1:(cut - 1)]))
set.seed(4711); R <- 40000L
cat("size of the CHI-SQUARED-calibrated test at nominal 0.05 (the way the paper uses B_n)\n")
cat("   n  stat     poly_p6 spline_p6   poly_p8 spline_p8\n")
for (n in c(100L, 200L)) {
  BB <- make_bases(n)
  M <- t(vapply(seq_len(R), function(i) all_stats(sample.int(n), sample.int(n), BB, n), numeric(20)))
  for (s in c("B", "P")) {
    f <- function(b, p) mean(M[, sprintf("%s_%s_%d", s, b, p)] > qchisq(0.95, p^2))
    cat(sprintf("%4d  %s_n      %.4f    %.4f      %.4f    %.4f\n", n, s,
                f("poly", 6), f("spline", 6), f("poly", 8), f("spline", 8)))
  }
}
cat("Monte Carlo s.e. about 0.0011\n")
