## =====================================================================
##  Calibration of the order scan.  Produces Tables S7 and S8 of the
##  Supplementary Material: the null correlation matrix of the Pillai
##  profile, and the threshold c_alpha against the Bonferroni threshold.
##
##  WHY TWO BATCHES.  The per-order null laws turn each P_n(k,k) into a
##  p-value; the law of min_k p_k / w_k fixes c_alpha.  Using ONE batch for
##  both would evaluate the scan on the very draws that defined its own
##  reference distribution, and the rule would over-reject.  The two
##  batches below are independent by construction.
##
##  Runs in about half a minute on one core.  Nothing here touches the
##  main power grid; it is a standalone calibration study.
## =====================================================================
library(dependence)
sp <- dependence:::scan_profile
sb <- dependence:::scan_basis

OUTDIR <- local({
  d <- tryCatch(dirname(rstudioapi::getSourceEditorContext()$path), error = function(e) NA)
  if (is.na(d) || !nzchar(d)) d <- getwd()
  p <- file.path(d, "results"); dir.create(p, showWarnings = FALSE); p
})

SEED <- 20261014L
NS   <- c(50, 100, 200)
KS   <- c(4, 5, 6)              # two orders above rule (12); see Section 2.6
R    <- 20000L                  # datasets per batch
ALPHA <- 0.05

set.seed(SEED)
rows <- list(); corrs <- list()
for (i in seq_along(NS)) {
  n <- NS[i]; K <- KS[i]; B <- sb(n, K, "poly")

  ## batch 1 -- the per-order null laws
  M1 <- matrix(0, R, K)
  for (r in seq_len(R)) M1[r, ] <- sp(runif(n), runif(n), K, "poly", "first", B)
  laws <- lapply(seq_len(K), function(k) sort(M1[, k]))

  ## batch 2 -- independent draws, converted to p-values
  M2 <- matrix(0, R, K)
  for (r in seq_len(R)) M2[r, ] <- sp(runif(n), runif(n), K, "poly", "first", B)
  P  <- sapply(seq_len(K), function(k) 1 - findInterval(M2[, k], laws[[k]]) / R)
  mp <- apply(P, 1, min)                      # uniform weights w_k = 1/K

  calpha <- as.numeric(quantile(mp, ALPHA, type = 1))
  rows[[i]] <- data.frame(n = n, K = K, c_alpha = calpha,
                          bonferroni = ALPHA / K,
                          ratio = calpha / (ALPHA / K),
                          size_bonferroni = mean(mp <= ALPHA / K))
  corrs[[i]] <- cor(M1)
  cat(sprintf("n=%3d K=%d  c_%.2f=%.5f  alpha/K=%.5f  ratio=%.2f  Bonferroni size=%.4f\n",
              n, K, ALPHA, calpha, ALPHA / K, calpha / (ALPHA / K),
              rows[[i]]$size_bonferroni))
}
tab <- do.call(rbind, rows)
names(corrs) <- paste0("n", NS)

write.csv(tab, file.path(OUTDIR, "scan_calibration.csv"), row.names = FALSE)
saveRDS(list(table = tab, corr = corrs),
        file.path(OUTDIR, "scan_calibration.rds"))
cat("\nnull correlations of the Pillai profile:\n")
for (i in seq_along(NS)) { cat("\nn =", NS[i], "\n"); print(round(corrs[[i]], 2)) }
message("Written: ", file.path(OUTDIR, "scan_calibration.csv"))
