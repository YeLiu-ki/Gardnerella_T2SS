#!/usr/bin/env Rscript

required_packages <- c("readr", "dplyr", "tidyr", "ggplot2", "patchwork", "ragg", "png")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Install required packages first: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(png)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) stop("Run this file with Rscript.")
script_dir <- dirname(normalizePath(sub("^--file=", "", script_arg)))
project_dir <- normalizePath(file.path(script_dir, ".."))
data_dir <- file.path(project_dir, "data")
output_dir <- file.path(project_dir, "output")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

bv_file <- file.path(data_dir, "BV_gardnerella_per_sample.tsv")
health_file <- file.path(data_dir, "Health_gardnerella_per_sample.tsv")
gp_filter_file <- file.path(data_dir, "02_filtered_samples_Gp.csv")
gv_filter_file <- file.path(data_dir, "02_filtered_samples_Gv.csv")
out_prefix <- file.path(output_dir, "gardnerella_abundance_filtered")

species_to_plot <- c("G_vaginalis", "G_piotii")
species_labels <- c(G_vaginalis = "G. vaginalis (Gv)", G_piotii = "G. piotii (Gp)")
species_colors <- c(G_vaginalis = "#0647C8", G_piotii = "#158B2C")
group_colors <- c(Health = "#0755D9", BV = "#FF3B3B")
group_order <- c("Health", "BV")

read_group <- function(path, expected_group) {
  x <- read_tsv(path, show_col_types = FALSE)
  required <- c("group", "sample", "species", "relative_abundance_pct")
  if (!all(required %in% names(x))) {
    stop("Missing columns in ", basename(path), ": ", paste(setdiff(required, names(x)), collapse = ", "))
  }
  if (!all(x$group == expected_group)) stop("Unexpected group in ", basename(path))
  x
}

abundance_raw <- bind_rows(
  read_group(health_file, "Health"),
  read_group(bv_file, "BV")
) %>%
  filter(species %in% species_to_plot) %>%
  transmute(
    group = factor(group, levels = group_order),
    sample,
    sample_key = sub("_run[0-9]+$", "", sample),
    species = factor(species, levels = species_to_plot),
    abundance = pmax(0, pmin(100, as.numeric(relative_abundance_pct)))
  )

read_filter <- function(path, species_name) {
  x <- read_csv(path, show_col_types = FALSE)
  required <- c("Sample", "Group")
  if (!all(required %in% names(x))) {
    stop("Missing columns in ", basename(path), ": ", paste(setdiff(required, names(x)), collapse = ", "))
  }
  x %>%
    transmute(species = factor(species_name, levels = species_to_plot),
              sample_key = as.character(Sample), filter_group = as.character(Group)) %>%
    distinct()
}

filter_map <- bind_rows(
  read_filter(gv_filter_file, "G_vaginalis"),
  read_filter(gp_filter_file, "G_piotii")
)
if (anyDuplicated(filter_map[c("species", "sample_key")])) stop("Duplicate filtered sample keys.")

dat <- abundance_raw %>% inner_join(filter_map, by = c("species", "sample_key"))
if (nrow(dat) != nrow(filter_map)) {
  missing_keys <- anti_join(filter_map, abundance_raw, by = c("species", "sample_key"))
  stop(nrow(missing_keys), " filtered species/sample entries were not found in abundance tables.")
}
if (any(as.character(dat$group) != dat$filter_group)) stop("Group labels disagree between inputs.")
dat <- dat %>% select(-filter_group)
if (any(count(dat, group, sample, species)$n != 1)) stop("Duplicate group/sample/species rows.")

get_n <- function(sp, grp) {
  dat %>% filter(species == sp, group == grp) %>%
    summarise(n = n_distinct(sample)) %>% pull(n)
}

theme_set(theme_classic(base_size = 9, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.4, colour = "black"),
        axis.ticks = element_line(linewidth = 0.4, colour = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5), panel.grid = element_blank()))

make_sample_panel <- function(sp, grp, show_x_title = FALSE) {
  wide <- dat %>%
    filter(species == sp, group == grp) %>%
    arrange(abundance, sample) %>%
    mutate(sample_index = row_number(), Others = 100 - abundance) %>%
    select(sample_index, abundance, Others) %>%
    pivot_longer(c(abundance, Others), names_to = "component", values_to = "value") %>%
    mutate(component = factor(component, levels = c("Others", "abundance")))

  ggplot(wide, aes(sample_index, value, fill = component)) +
    geom_col(width = 0.94, colour = NA) +
    scale_fill_manual(values = c(abundance = unname(species_colors[sp]), Others = "#D9D9D9"),
                      breaks = c("abundance", "Others"),
                      labels = c(species_labels[[sp]], "Others")) +
    scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100), expand = c(0, 0)) +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = sprintf("%s (n = %d)", ifelse(grp == "Health", "Healthy", "BV"), get_n(sp, grp)),
         x = if (show_x_title) "Samples (sorted by relative abundance)" else NULL,
         y = if (grp == "Health") "Relative abundance (%)" else NULL, fill = NULL) +
    theme(plot.title = element_text(size = 10, colour = unname(group_colors[grp]), margin = margin(3, 3, 3, 3)),
          plot.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(colour = "#BEBEBE", fill = NA, linewidth = 0.3),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_text(size = 9, face = "bold"), axis.title.y = element_text(size = 9),
          legend.position = if (grp == "BV") "right" else "none",
          legend.text = element_text(size = 8), legend.key.height = grid::unit(4, "mm"))
}

p_to_stars <- function(p) {
  if (is.na(p)) "NA" else if (p < 0.0001) "****" else if (p < 0.001) "***" else if (p < 0.01) "**" else if (p < 0.05) "*" else "ns"
}

make_summary_panel <- function(sp) {
  d <- dat %>% filter(species == sp)
  stats <- d %>% group_by(group) %>%
    summarise(n = n(), mean = mean(abundance), se = sd(abundance) / sqrt(n),
              ci = qt(0.975, df = pmax(n - 1, 1)) * se, .groups = "drop")
  wt <- wilcox.test(abundance ~ group, data = d, exact = FALSE)
  ymax_data <- max(d$abundance, stats$mean + stats$ci, na.rm = TRUE)
  bracket_y <- ymax_data * 1.12 + 0.5
  label_y <- bracket_y + ymax_data * 0.06 + 0.15
  upper <- label_y + ymax_data * 0.12 + 0.5

  ggplot(d, aes(group, abundance, colour = group)) +
    geom_col(data = stats, aes(group, mean, fill = group), width = 0.55, alpha = 0.28,
             colour = "black", linewidth = 0.45, inherit.aes = FALSE) +
    geom_errorbar(data = stats, aes(group, ymin = pmax(0, mean - ci), ymax = mean + ci),
                  width = 0.24, linewidth = 0.65, colour = "black", inherit.aes = FALSE) +
    geom_point(position = position_jitter(width = 0.10, height = 0, seed = 260622), size = 0.8, alpha = 0.62) +
    annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y, colour = "black", linewidth = 0.4) +
    annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - 0.04 * upper, colour = "black", linewidth = 0.4) +
    annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - 0.04 * upper, colour = "black", linewidth = 0.4) +
    annotate("text", x = 1.5, y = label_y, label = p_to_stars(wt$p.value), size = 3.0, fontface = "bold") +
    geom_text(data = stats, aes(group, y = mean + ci + 0.06 * upper,
                                label = sprintf("%.2f%%", mean), colour = group),
              size = 3.0, fontface = "bold", inherit.aes = FALSE) +
    scale_colour_manual(values = group_colors) + scale_fill_manual(values = group_colors) +
    scale_x_discrete(labels = c(Health = sprintf("Healthy\n(n = %d)", get_n(sp, "Health")),
                                BV = sprintf("BV\n(n = %d)", get_n(sp, "BV")))) +
    scale_y_continuous(limits = c(0, upper), expand = expansion(mult = c(0, 0.01))) +
    labs(title = species_labels[[sp]], x = NULL, y = "Relative abundance (%)") +
    theme(plot.title = element_text(size = 11, face = "bold.italic", colour = species_colors[[sp]], margin = margin(b = 6)),
          axis.text.x = element_text(size = 9, colour = "black"), axis.title.y = element_text(size = 9), legend.position = "none")
}

make_species_row <- function(sp, show_x_title = FALSE) {
  left <- (make_sample_panel(sp, "Health", show_x_title) | make_sample_panel(sp, "BV", show_x_title)) +
    plot_layout(widths = c(1, 1.32)) +
    plot_annotation(title = species_labels[[sp]],
                    theme = theme(plot.title = element_text(size = 12, face = "bold.italic", hjust = 0.5,
                                                            colour = species_colors[[sp]], margin = margin(b = 1))))
  (wrap_elements(full = left) | make_summary_panel(sp)) + plot_layout(widths = c(3.45, 1.2))
}

fig <- wrap_elements(full = make_species_row("G_vaginalis", FALSE)) /
  wrap_elements(full = make_species_row("G_piotii", TRUE)) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Per-sample and average relative abundance",
                  subtitle = "Bars: mean +/- 95% CI; points: individual samples; Wilcoxon rank-sum test",
                  theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                                plot.subtitle = element_text(size = 8.5, hjust = 0.5, colour = "#444444")))

ragg::agg_png(paste0(out_prefix, ".png"), width = 11.0, height = 7.25,
              units = "in", res = 600, background = "white")
print(fig)
dev.off()

# Embed the verified 600-dpi rendering in PDF. This avoids device-specific
# nested-patchwork rendering differences across operating systems.
rendered_figure <- png::readPNG(paste0(out_prefix, ".png"))
grDevices::cairo_pdf(paste0(out_prefix, ".pdf"), width = 11.0, height = 7.25)
grid::grid.newpage()
grid::grid.raster(rendered_figure, width = grid::unit(1, "npc"),
                  height = grid::unit(1, "npc"), interpolate = TRUE)
dev.off()

summary_table <- dat %>% group_by(species, group) %>%
  summarise(n = n(), mean_pct = mean(abundance),
            ci95_lower_pct = mean(abundance) - qt(0.975, df = n() - 1) * sd(abundance) / sqrt(n()),
            ci95_upper_pct = mean(abundance) + qt(0.975, df = n() - 1) * sd(abundance) / sqrt(n()),
            .groups = "drop")
write_tsv(summary_table, paste0(out_prefix, "_summary.tsv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo_gardnerella.txt"))
message("Written: ", paste0(out_prefix, c(".pdf", ".png", "_summary.tsv"), collapse = "\n         "))
