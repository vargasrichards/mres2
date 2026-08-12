lloq_pcr   <- 3      # log10 copies/mL (viral RNA / PCR)
lloq_titre <- 1.27   # log10 FFU/mL (infectious virus titre)
llod       <- 0      # log10, common to both assays
lloq_ifnalpha <- -0.31

simulx_seed  <- 20
monolix_seed <- simulx_seed

suppressPackageStartupMessages(library(lixoftConnectors))
suppressPackageStartupMessages(library(RsSimulx))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(forcats))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(here))

here::i_am("code/make_early_figs.R")

source(here::here("code", "paths.R"))
source(here::here("code", "simulation_utils.R"))

required_objects <- c("outputs_wanted", "lloq_pcr", "lloq_titre",
                      "lloq_ifnalpha", "llod", "remap_ids", "write_outputs",
                      "longify", "logtransform_select", "plot_vars")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1))]
if (length(missing_objects)) {
  stop("simulation_utils.R did not define: ",
       paste(missing_objects, collapse = ", "))
}

runest <- FALSE

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

stream_cols <- c(
  "PCR viral load"   = unname(okabe_ito_pal["blue"]),
  "Infectious titre" = unname(okabe_ito_pal["bluish_green"]),
  "IFN-α2a conc." = unname(okabe_ito_pal["vermillion"])
)

col_pcr   <- stream_cols[["PCR viral load"]]
col_titre <- stream_cols[["Infectious titre"]]
col_ifn   <- stream_cols[["IFN-α2a conc."]]

aux_cols <- c(
  "Infected cells"           = unname(okabe_ito_pal["orange"]),
  "IFN"                      = col_ifn,
  "log(Infected/interferon)" = unname(okabe_ito_pal["purple"])
)

state_cols <- c(
  T = unname(okabe_ito_pal["sky_blue"]),
  E = unname(okabe_ito_pal["yellow"]),
  I = unname(okabe_ito_pal["orange"]),
  R = unname(okabe_ito_pal["bluish_green"]),
  D = unname(okabe_ito_pal["grey"])
)

id_cols <- unname(okabe_ito_pal[c("blue", "vermillion", "bluish_green",
                                  "purple", "orange", "sky_blue")])

alpha_model       <- 1
alpha_cens_hlines <- 0.6
pointsiz          <- 0.9

theme_set(theme_bw(base_size = 14))

initializeLixoftConnectors(software = "monolix", force = TRUE)
loadProject(path_mlxtran_3_explicit)

original_data_filename <- getData()$dataFile
obsInfo <- getObservationInformation()

if (runest) {
  setPopulationParameterEstimationSettings(seed = monolix_seed)
  runPopulationParameterEstimation()
  runConditionalDistributionSampling()
  runConditionalModeEstimation()
  runStandardErrorEstimation()
  runLogLikelihoodEstimation()
}

indivParams <- getEstimatedIndividualParameters()
ind <- indivParams$conditionalMode

initializeLixoftConnectors(software = "simulx", force = TRUE)

if (runest || !file.exists(path_indi_remapped)) {
  res_indiv <- RsSimulx::simulx(
    project         = path_mlxtran_3_explicit,
    parameter       = ind,
    output          = outputs_wanted,
    group           = list(size = nrow(ind), level = "individual"),
    settings        = list(seed = simulx_seed, exportData = TRUE),
    saveSmlxProject = path_smlx_indiv
  )
  write_outputs(res_indiv, path_outputs_indiv)
  res_indi_remapped <- lapply(res_indiv, remap_ids)
  saveRDS(res_indi_remapped, path_indi_remapped)
} else {
  res_indi_remapped <- readRDS(path_indi_remapped)
}

sim_ids <- res_indi_remapped$ValueVt$id
id_levels <- if (is.factor(sim_ids)) levels(sim_ids) else sort(unique(as.character(sim_ids)))

harmonise_ids <- function(df, levels_vec) {
  df$id <- factor(as.character(df$id), levels = levels_vec)
  df
}

obs_y1 <- harmonise_ids(obsInfo$y1, id_levels)
obs_y2 <- harmonise_ids(obsInfo$y2, id_levels)
obs_y3 <- harmonise_ids(obsInfo$y3, id_levels) |>
  dplyr::mutate(y3_censored = y3 < lloq_ifnalpha,
                y3 = ifelse(y3 < lloq_ifnalpha, lloq_ifnalpha, y3))

unmatched <- sum(is.na(obs_y1$id))
if (unmatched > 0) {
  stop(unmatched, " observation rows have ids not present in the simulated ",
       "output. Check that remap_ids() and the Monolix ids agree before plotting.")
}

obs_vi_vt <- dplyr::left_join(obs_y1, obs_y2,
                              by = c("id", "time"))

indiv_both <- dplyr::left_join(
  res_indi_remapped$ValueVi,
  res_indi_remapped$ValueVt,
  by = c("id", "time")
)

if (!exists("extended")) extended <- FALSE

if (extended) {

p_dual <- ggplot(indiv_both) +
  geom_line(aes(x = ValueVt, y = ValueVi), colour = col_pcr) +
  geom_point(data = obs_vi_vt, aes(x = y1, y = y2), size = pointsiz) +
  facet_wrap(~id, scales = "fixed") +
  labs(y = "Infectious titre (log10 FFU/mL)",
       x = "Viral load (log10 copies/mL)",
       title = "Relation between viral load and infectious virus")
save_fig(p_dual, "dual_plot.pdf")

ppcr <- ggplot() +
  geom_line(data = res_indi_remapped$ValueVt, aes(time, ValueVt),
            colour = col_pcr) +
  geom_point(data = obs_y1, aes(time, y1), size = .8) +
  facet_wrap(~id, scales = "fixed") +
  labs(x = "Time /days post-inoculation", y = "log10(copies/mL)",
       title = "Viral load: fits of cellular immunity with explicit IFN")
save_fig(ppcr, "pcr_fits_cell_immun.pdf")

pffu <- ggplot() +
  geom_line(data = res_indi_remapped$ValueVi, aes(time, ValueVi),
            colour = col_titre) +
  geom_point(data = obs_y2, aes(time, y2), size = .8) +
  facet_wrap(~id, scales = "fixed") +
  labs(x = "Time /days post-inoculation", y = "log10(FFU/mL)",
       title = "Infectious titre: fits of cellular immunity with explicit IFN")
save_fig(pffu, "ffu_fits_cell_immun.pdf")

pifn <- ggplot() +
  geom_line(data = res_indi_remapped$ValueF, aes(time, ValueF),
            colour = col_ifn) +
  geom_point(data = obs_y3, aes(time, y3), size = .8) +
  facet_wrap(~id, scales = "fixed") +
  labs(x = "Time /days post-inoculation", y = "log10 IFN-α2a (pg/mL)",
       title = "Interferon: fits of cellular immunity with explicit IFN")
save_fig(pifn, "ifn_fits_cell_immun.pdf")

}

thresholds <- tibble::tibble(
  yint = c(lloq_pcr, lloq_titre, lloq_ifnalpha, llod),
  label = factor(
    c("LLOQ PCR (3 log10 copies/mL)",
      "LLOQ titre (1.27 log10 FFU/mL)",
      "LLOQ IFN-α2a (-0.3 log10 pg/mL)",
      "LLOD (0 log10, both assays)"),
    levels = c("LLOQ PCR (3 log10 copies/mL)",
               "LLOQ titre (1.27 log10 FFU/mL)",
               "LLOQ IFN-α2a (-0.3 log10 pg/mL)",
               "LLOD (0 log10, both assays)")
  )
)

pcr_string <- "Viral load (PCR)"

p_tristream <- ggplot() +
  geom_hline(
    data = thresholds,
    aes(yintercept = yint, linetype = label),
    colour = "grey35", linewidth = 0.35,
    alpha = alpha_cens_hlines, inherit.aes = FALSE
  ) +
  geom_line(data = res_indi_remapped$ValueVt,
            aes(time, ValueVt, colour = pcr_string), alpha = alpha_model) +
  geom_line(data = res_indi_remapped$ValueVi,
            aes(time, ValueVi, colour = "Infectious titre"), alpha = alpha_model) +
  geom_line(data = res_indi_remapped$ValueF,
            aes(time, ValueF, colour = "IFN-α2a conc."), alpha = alpha_model) +
  geom_point(data = obs_y1,
             aes(time, y1, colour = pcr_string,
                 shape = y1 <= lloq_pcr), size = pointsiz) +
  geom_point(data = obs_y2,
             aes(time, y2, colour = "Infectious titre",
                 shape = y2 <= lloq_titre), size = pointsiz) +
  geom_point(data = obs_y3,
             aes(time, y3, colour = "IFN-α2a conc.",
                 shape = y3 <= lloq_ifnalpha), size = pointsiz) +
  scale_colour_manual(
    values = stream_cols,
    breaks = names(stream_cols)
  ) +
  scale_linetype_manual(
    values = c("LLOQ PCR (3 log10 copies/mL)"        = "dashed",
               "LLOQ titre (1.27 log10 FFU/mL)"      = "dotted",
               "LLOQ IFN-α2a (-0.3 log10 pg/mL)" = "dotdash",
               "LLOD (0 log10, both assays)"         = "longdash"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 1),
    labels = c("FALSE" = "Quantifiable", "TRUE" = "Censored (≤ LLOQ)"),
    name = NULL,
    guide = "none"
  ) +
  facet_wrap(~id, scales = "fixed") +
  labs(
    x = "Time /days post-inoculation",
    y = "log10 concentration (/mL)",
    colour = NULL,
    title = "Viral load, infectious titre, and IFN-α2a: data and model fit"
  ) +
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(
    colour = guide_legend(order = 1, override.aes = list(shape = NA))
  )

save_fig(p_tristream, "Fig_2_pcr_ffu_ifn_fits_cell_immun.pdf",
         width = 12, height = 10)

res_indi_remapped$T_freq <- res_indi_remapped$T_freq |>
  dplyr::mutate(Rfreq = 1 - T_freq)

frac_refrac <- ggplot(res_indi_remapped$T_freq, aes(time, Rfreq)) +
  geom_line(colour = unname(okabe_ito_pal["bluish_green"])) +
  facet_wrap(~id, scales = "fixed") +
  labs(x = "Time /days post-inoculation",
       y = "Fraction of target cells refractory",
       title = "Refractory fraction")
save_fig(frac_refrac, "frac_refrac.pdf")

res_indi_remapped$I <- res_indi_remapped$I |> dplyr::mutate(logI = log10(I + 1))

res_indi_remapped$logI_minus_logF <- dplyr::left_join(
  res_indi_remapped$ValueF,
  res_indi_remapped$I,
  by = c("id", "time")
) |>
  dplyr::mutate(logI_minus_logF = logI - ValueF) |>
  dplyr::select(!c(logI, I, ValueF))

combined_indiv <- dplyr::bind_rows(lapply(res_indi_remapped, longify)) |>
  logtransform_select()

target_cells <- plot_vars(combined_indiv, "T", ystring = "log10 Target cells (T)") +
  scale_colour_manual(values = c("T" = state_cols[["T"]]))
save_fig(target_cells, "target_cells.pdf")

combined_indiv_ii <- combined_indiv |>
  dplyr::filter(variable %in% c("logI", "ValueF", "logI_minus_logF")) |>
  dplyr::mutate(variable = dplyr::case_when(
    variable == "logI"             ~ "Infected cells",
    variable == "ValueF"           ~ "IFN",
    variable == "logI_minus_logF"  ~ "Infected/interferon"
  ))

infected_and_interferon <- plot_vars(
  combined_indiv_ii,
  c("Infected cells", "IFN", "Infected/interferon"),
  ystring = "log10 concentration (/mL)",
  plot_title = "Interferon and infected cell dynamics: individual fits"
) +
  geom_hline(yintercept = lloq_ifnalpha, linetype = "dotted",
             colour = "grey50", linewidth = 0.5) +
  geom_point(
    data = obs_y3,
    aes(x = time, y = y3),
    colour = col_ifn, inherit.aes = FALSE, size = pointsiz
  ) + scale_x_continuous(limits = c(0, 15)) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", legend.title = element_blank())

save_fig(infected_and_interferon, "Fig_5_infected_vs_ifn.pdf")

infected_cells <- ggplot(res_indi_remapped$I, aes(time, logI)) +
  geom_line(colour = aux_cols[["Infected cells"]]) +
  facet_wrap(~id, scales = "fixed") +
  labs(x = "Time /days post-inoculation",
       y = "log10 infected cell concentration")
save_fig(infected_cells, "I.pdf")

p_T_to_R <- ggplot(res_indi_remapped$flow_T_to_R) +
  geom_line(aes(time, flow_T_to_R), colour = unname(okabe_ito_pal["purple"])) +
  facet_wrap(~id, scales = "fixed") +
  labs(x = "Time /days post-inoculation",
       y = "Flow T → R",
       title = "Target to refractory cell flow rate")
save_fig(p_T_to_R, "flow_T_to_R.pdf")

rebound_ids <- c("13", "15")

finite_diff <- function(df, ids, var, new_name) {
  df |>
    dplyr::filter(id %in% ids, variable == var) |>
    dplyr::arrange(id, time) |>
    dplyr::group_by(id) |>
    dplyr::mutate(value = c(NA_real_, diff(value) / diff(time)),
                  variable = new_name) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(value))
}

ddt_vt <- finite_diff(combined_indiv, rebound_ids, "ValueVt", "ddt_Vt")
ddt_f  <- finite_diff(combined_indiv, rebound_ids, "ValueF",  "ddt_F")

rebound_windows <- ddt_vt |>
  dplyr::arrange(id, time) |>
  dplyr::group_by(id) |>
  dplyr::reframe({
    t <- time; s <- sign(value)
    if (length(s) < 2) {
      tibble::tibble(start = NA_real_, end = NA_real_)
    } else {
      pos_to_neg <- which(s[-length(s)] > 0 & s[-1] < 0) + 1
      neg_to_pos <- which(s[-length(s)] < 0 & s[-1] > 0) + 1

      peak_idx  <- if (length(pos_to_neg)) pos_to_neg[1] else NA_integer_
      start_idx <- if (!is.na(peak_idx)) {
        cand <- neg_to_pos[neg_to_pos > peak_idx]
        if (length(cand)) cand[1] else NA_integer_
      } else NA_integer_
      end_idx <- if (!is.na(start_idx)) {
        cand <- pos_to_neg[pos_to_neg > start_idx]
        if (length(cand)) cand[1] else NA_integer_
      } else NA_integer_

      tibble::tibble(
        start = if (!is.na(start_idx)) t[start_idx] else NA_real_,
        end   = if (!is.na(end_idx)) t[end_idx]
        else if (!is.na(start_idx)) max(t) else NA_real_
      )
    }
  }) |>
  tidyr::drop_na(start, end) |>
  dplyr::mutate(id = factor(id, levels = rebound_ids))

if (nrow(rebound_windows) == 0) {
  warning("No rebound window detected for ids: ",
          paste(rebound_ids, collapse = ", "))
}

first_crossing <- function(df, direction = c("neg_to_pos", "pos_to_neg"),
                           t_min = -Inf) {
  direction <- match.arg(direction)
  df <- df |> dplyr::filter(time > t_min) |> dplyr::arrange(time)
  if (nrow(df) < 2) return(NA_real_)
  s <- sign(df$value)
  hits <- if (direction == "neg_to_pos") {
    which(s[-length(s)] < 0 & s[-1] > 0)
  } else {
    which(s[-length(s)] > 0 & s[-1] < 0)
  }
  if (!length(hits)) return(NA_real_)
  df$time[hits[1] + 1]
}

event_markers <- purrr::map_dfr(rebound_ids, function(this_id) {
  ct <- combined_indiv |>
    dplyr::filter(id == this_id, variable == "changeinT") |>
    dplyr::arrange(time)
  ft <- ddt_f |>
    dplyr::filter(id == this_id)

  tibble::tibble(
    id    = this_id,
    event = c("target_cells_increase", "ifn_decline"),
    time  = c(first_crossing(ct, "neg_to_pos"),
              first_crossing(ft, "pos_to_neg"))
  )
}) |>
  dplyr::mutate(id = factor(id, levels = rebound_ids)) |>
  tidyr::drop_na(time)

plot_df <- combined_indiv |>
  dplyr::filter(id %in% rebound_ids,
                variable %in% c("changeinT", "ValueF", "ValueVt")) |>
  dplyr::bind_rows(ddt_vt) |>
  dplyr::mutate(
    variable = factor(variable,
                      levels = c("ddt_Vt", "ValueVt", "changeinT", "ValueF"),
                      labels = c("d/dt log10 Vt", "log10 Vt", "Δ target cells",
                                 "log10 IFN-α2a")),
    id = factor(as.character(id), levels = rebound_ids)
  )

rebound_analysis <- ggplot(plot_df, aes(time, value, colour = id)) +
  geom_rect(
    data = rebound_windows,
    aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
    fill = "grey88", inherit.aes = FALSE
  ) +
  geom_line(linewidth = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.3) +
  geom_vline(
    data = event_markers,
    aes(xintercept = time, linetype = event),
    colour = "turquoise4", linewidth = 0.9, inherit.aes = FALSE
  ) +
  scale_linetype_manual(
    values = c(target_cells_increase = "dotted", ifn_decline = "dotdash"),
    labels = c(target_cells_increase = "Target cells begin increasing",
               ifn_decline = "IFN begins decreasing"),
    name = NULL
  ) + scale_colour_manual(values = id_cols, name = "ID") +
  facet_grid(rows = vars(variable), cols = vars(id), scales = "free_y") +
  labs(x = "Time /days post-inoculation", y = "Concentration /mL") +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom", axis.text = element_text(size = 14))

save_fig(rebound_analysis, "Fig_4_rebound.pdf", width = 10, height = 11)
