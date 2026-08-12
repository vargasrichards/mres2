outputs_wanted <- list(
  name = c("ValueVt",
           "ValueVi",
           "ValueF",
           "flow_T_to_R",
           "flow_R_to_T",
           "T_freq",
           "changeinT",
           "T","D",
           "I", "foi", "E", "R", "delta"),
  time = seq(0, 20, by = 0.1)
)

is_output_df <- function(df) {
  if (!all(c("id", "time") %in% names(df))) return(FALSE)
  value_cols <- setdiff(names(df), c("id", "time"))
  length(value_cols) > 0 && all(vapply(df[value_cols], is.numeric, logical(1)))
}

longify <- function(df) {
  if (!is_output_df(df)) return(NULL)
  tidyr::pivot_longer(df, cols = -c(id, time), names_to = "variable", values_to = "value")
}

write_outputs <- function(res_list, outdir) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  for (nm in names(res_list)) {
    write.csv(res_list[[nm]], file = file.path(outdir, paste0(nm, ".csv")), row.names = FALSE)
  }
}

remap_ids <- function(df) {
  if ("original_id" %in% names(df)) {
    df |>
      dplyr::select(-id) |>
      dplyr::rename(id = original_id) |>
      dplyr::mutate(
        id = factor(as.character(id),
                    levels = as.character(
                      sort(as.numeric(
                        as.character(unique(id))))))
      )
  } else {
    df
  }
}

plot_vars <- function(df,
                      vars,
                      plot_title = NULL,
                      ystring = NULL,
                      xstring = "Time post-challenge /days",
                      lower_y_lim = FALSE) {
  p <- df |>
    dplyr::filter(variable %in% vars) |>
    ggplot(aes(time, value, colour = variable, group = interaction(id, variable))) +
    geom_line() +
    facet_wrap(~id, scales = "fixed") + labs(title = plot_title, x = xstring, y = ystring) +
    theme_bw() + theme(text = element_text(size = 14), strip.text.x = element_text(size = 14),
                            axis.text = element_text(size = 14))

  if (lower_y_lim != FALSE) {
    p <- p + lims(y = c(0, NA))
  }
  p
}

logtransform_select <- function(df, vars_to_transform = c("T", "flow_R_to_T", "flow_T_to_R", "foi")) {
    dft <- df |>
  dplyr::mutate(value = dplyr::if_else(variable %in% vars_to_transform,
                                           log10(value),
                                           value))

  dft

}

