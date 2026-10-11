library(tidyverse)

res <- data.table::fread("inst/paper/results/power_continuous.csv")

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

ggsave("inst/paper/results/power_continuous.pdf", width = 6.5, height = 8.5, device = cairo_pdf)


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
