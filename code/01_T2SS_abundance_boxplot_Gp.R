#!/usr/bin/env Rscript

required_packages <- c("ggplot2", "ggpubr", "dplyr", "readr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(readr)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
project_dir <- normalizePath(file.path(script_dir, ".."))
input_file <- file.path(project_dir, "data", "02_filtered_samples_Gp.csv")
output_dir <- file.path(project_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    Group = if_else(Group == "Healthy", "Health", Group),
    Group = factor(Group, levels = c("Health", "BV"))
  )

required_columns <- c("Sample", "Completeness", "Group", "RPM")
if (!all(required_columns %in% names(dat))) {
  stop("Input is missing: ", paste(setdiff(required_columns, names(dat)), collapse = ", "))
}
if (anyNA(dat[required_columns])) stop("Required columns contain missing values.")
if (any(dat$RPM <= 0)) stop("RPM must be positive before log10 transformation.")
if (any(dat$Completeness <= 0.6)) warning("Input contains rows with Completeness <= 0.6.")

set.seed(260622)
group_colors <- c(BV = "#D55E00", Health = "#0072B2")
p <- ggplot(dat, aes(x = Group, y = log10(RPM), fill = Group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.75,
               color = "black", linewidth = 0.6) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.5) +
  scale_fill_manual(values = group_colors) +
  theme_classic(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(linewidth = 0.8),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  ) +
  labs(y = "Log 10 (T2SS abundance (RPM))", x = NULL) +
  stat_compare_means(
    comparisons = list(c("BV", "Health")),
    method = "wilcox.test",
    label = "p.format"
  )

pdf_file <- file.path(output_dir, "T2SS_abundance_boxplot_Gp.pdf")
ggsave(pdf_file, p, width = 3, height = 4, units = "in", device = cairo_pdf)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_Gp.txt"))
message("Written: ", pdf_file)
