library(here)

dir_data     <- here::here("data")
dir_models   <- here::here("models")
dir_struct   <- here::here("models", "structural_models")
dir_figures  <- here::here("figures")
dir_outputs  <- here::here("outputs")

path_data            <- file.path(dir_data, "chim_filtered_3stream.csv")
path_mlxtran_3_explicit <- file.path(dir_models, "3_explicit.mlxtran")
path_mlxtran_3_prop     <- file.path(dir_models, "3_prop.mlxtran")
path_mlxtran_2_explicit <- file.path(dir_models, "2_explicit.mlxtran")
path_mlxtran_2_prop     <- file.path(dir_models, "2_prop.mlxtran")
path_struct_explicit <- file.path(dir_struct, "explicit_IFN.txt")
path_struct_prop     <- file.path(dir_struct, "prop_IFN.txt")

path_smlx_indiv      <- file.path(dir_outputs, "simulx", "individual_fits.smlx")
path_outputs_indiv   <- file.path(dir_outputs, "simulx_outputs_indiv")
path_indi_remapped   <- file.path(dir_outputs, "indi_remapped.RDS")

for (d in c(dir_figures, dir_outputs,
            dirname(path_smlx_indiv), path_outputs_indiv)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

fig_path <- function(filename) file.path(dir_figures, filename)

save_fig <- function(plot, filename, width = 12, height = 9, ...) {
  ggplot2::ggsave(
    filename = fig_path(filename), plot = plot,
    width = width, height = height, device = cairo_pdf, ...
  )
}
