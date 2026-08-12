library(ggplot2)
library(dplyr)
library(tidyr)
library(here)

here::i_am("code/Fig_3_cell_state_plot.R")

source(here::here("code", "paths.R"))
source(here::here("code", "simulation_utils.R"))

okabe_ito_pal <- c(
  black        = "#000000",
  orange       = "#E69F00",
  sky_blue     = "#56B4E9",
  bluish_green = "#009E73",
  yellow       = "#F0E442",
  blue         = "#0072B2",
  vermillion   = "#D55E00",
  purple       = "#CC79A7",
  grey         = "#999999"
)

state_cols <- c(
  T = unname(okabe_ito_pal["sky_blue"]),
  E = unname(okabe_ito_pal["yellow"]),
  I = unname(okabe_ito_pal["orange"]),
  R = unname(okabe_ito_pal["bluish_green"]),
  D = unname(okabe_ito_pal["grey"])
)

theme_set(theme_bw(base_size = 14))

if (!file.exists(path_indi_remapped)) {
  stop("Simulated output not found at ", path_indi_remapped,
       ". Run code/make_early_figs.R (or run_all.R) first.")
}

res_indi_remapped <- readRDS(path_indi_remapped)

combined_indiv_raw <- dplyr::bind_rows(lapply(res_indi_remapped, longify))

compartments_wide <- combined_indiv_raw |>
  dplyr::filter(variable %in% c("D", "E", "I", "R", "T")) |>
  tidyr::pivot_wider(
    id_cols     = c(id, time),
    names_from  = variable,
    values_from = value,
    values_fn   = mean
  )

compartments_wide <- compartments_wide |>
  dplyr::mutate(N = D + E + I + R + T,
                EI_prop = (E + I) / N)

compartments_long <- compartments_wide |>
  dplyr::mutate(dplyr::across(c(D, E, I, R, T), ~ .x / N)) |>
  dplyr::select(id, time, D, E, I, R, T) |>
  tidyr::pivot_longer(cols = c(D, E, I, R, T),
                      names_to = "state", values_to = "value_n") |>
  dplyr::mutate(state = factor(state, levels = c("T", "E", "I", "R", "D")))

max_EI_prop <- max(compartments_wide$EI_prop, na.rm = TRUE)
scale_factor <- if (max_EI_prop > 0) 1 / max_EI_prop else 1

state_analysis <- ggplot() +
  geom_area(data = compartments_long, aes(time, value_n, fill = state)) +
  geom_line(data = compartments_wide, aes(time, EI_prop * scale_factor),
            colour = "black", linetype = "dashed", linewidth = 0.5) +
  facet_wrap(vars(id)) +
  scale_fill_manual(values = state_cols, breaks = c("T", "E", "I", "R", "D")) +
  scale_y_continuous(
    name = "Proportion of cells in state",
    sec.axis = sec_axis(~ . / scale_factor, name = "Proportion E + I (dashed)")
  ) +
  labs(x = "Time /days post-inoculation", fill = "State") +
  theme(
    text = element_text(size = 18),
    strip.text = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 16),
    legend.title = element_text(size = 18, face = "bold")
  )

save_fig(state_analysis, "Fig_3_cell_state_analysis.pdf")
