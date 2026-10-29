#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this file with Rscript.")
project_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))

scripts <- file.path(project_dir, "code", c(
  "01_T2SS_abundance_boxplot_Gp.R",
  "02_T2SS_abundance_boxplot_Gv.R",
  "03_gardnerella_abundance_filtered.R"
))

for (script in scripts) {
  message("Running ", basename(script), " ...")
  status <- system2(file.path(R.home("bin"), "Rscript"), script)
  if (status != 0) stop("Failed: ", basename(script))
}
message("All figures were reproduced in: ", file.path(project_dir, "output"))
