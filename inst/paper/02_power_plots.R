library(tidyverse)

## ---------------------------- where output goes -----------------------
## EVERYTHING this script writes -- cached cells, critical-value tables,
## the assembled CSV, the LaTeX fragments, the figures -- goes to
## inst/paper/results inside the package, whatever R's working directory
## happens to be.
##
## That was not true before, and it is why results appeared in two places.
## OUTDIR was the relative string "results", so a run started at the
## package root wrote cells to dependence/results while a run started in
## inst/paper wrote them to inst/paper/results; meanwhile OUTSTEM carried no
## directory at all, so the assembled CSV and the LaTeX fragments went to
## the working directory rather than into OUTDIR, and the cells and the
## tables built from them could end up in different folders.  The three
## companion scripts each had a third convention.  All of them now use this
## block, so there is exactly one results directory.
##
## The directory is found from the location of THIS FILE, which Rscript,
## source() and the RStudio editor each report differently, with the
## working directory as the last resort -- that is the case of "Run
## selection", where nothing identifies the file.  From each candidate the
## search walks up to the package root (the DESCRIPTION file) and back down
## to inst/paper, and also accepts a directory that IS the paper directory,
## which is what an installed copy of the package looks like.
paper_dir <- function() {
  cand <- character(0)
  a <- commandArgs(trailingOnly = FALSE)                      # Rscript
  m <- grep("^--file=", a, value = TRUE)
  if (length(m)) cand <- c(cand, dirname(sub("^--file=", "", m[1L])))
  for (i in rev(seq_len(sys.nframe()))) {                     # source()
    of <- tryCatch(sys.frame(i)$ofile, error = function(e) NULL)
    if (is.character(of) && length(of) == 1L) {
      cand <- c(cand, dirname(of)); break
    }
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) &&       # RStudio editor
      isTRUE(tryCatch(rstudioapi::isAvailable(), error = function(e) FALSE))) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (is.character(p) && length(p) == 1L && nzchar(p)) cand <- c(cand, dirname(p))
  }
  cand <- c(cand, getwd())                                    # last resort
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
  stop("cannot locate inst/paper from \"", getwd(), "\". ",
       "Set OUTDIR to the results directory by hand just below this block.",
       call. = FALSE)
}

PAPER   <- paper_dir()
OUTDIR  <- file.path(PAPER, "results")
OUTSTEM <- file.path(OUTDIR, "power_continuous")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
message("results directory: ", OUTDIR)


res <- data.table::fread(paste0(OUTSTEM, ".csv"))

# PLOTS

res |>
  filter(q != 1, test != "P_n", scheme != "linear") |>
  rename(`sd(noise)` = nu) |>
  mutate(scheme = recode(scheme,
                         "heteroscedastic" = "heterosc.",
                         "fourth root" = "4th root",
                         "two curves" = "2 curves"),
         n = recode(as.factor(n),
                    "50"  = "n = 50",
                    "100" = "n = 100",
                    "200" = "n = 200")) |>
  ggplot(aes(`sd(noise)`, y = power, color = test, linetype = test)) +
  geom_line(linewidth = .8) +
  facet_grid(scheme~n) +
  theme_bw()

ggsave(paste0(OUTSTEM, ".pdf"), width = 6.5, height = 8.5, device = cairo_pdf)


# q = 1

# TABLES


res |>
  filter(q != 1, scheme != "linear") |>
  rename(`sd(noise)` = nu) |>
  mutate(scheme = recode(scheme,
                         "heteroscedastic" = "heterosc.",
                         "fourth root" = "4th root",
                         "two curves" = "2 curves"),
         n = recode(as.factor(n),
                    "50"  = "n = 50",
                    "100" = "n = 100",
                    "200" = "n = 200")) |>
  group_by(n, test) |>
  summarise(mean_power = mean(power)) |>
  ungroup() |>
  pivot_wider(names_from = test, values_from = mean_power) |>
  xtable::xtable()
