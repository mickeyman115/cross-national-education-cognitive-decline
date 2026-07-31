###############################################################################
# Submission-ready manuscript figures from locked aggregate outputs
#
# Figure contract
# Core conclusion: the intermediate-versus-low education difference in 0-4-year
# memory change is larger in CHARLS than in HRS, ELSA, or supported SHARE
# countries, and this comparison is stable to measured first-return selection.
# Archetype: quantitative comparison series with one study-flow schematic.
# Export: 183-mm figures; editable SVG/PDF and 600-dpi TIFF.
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

source("analysis/00_config.R")
proj_dir <- PROJECT_DIR
out_dir <- file.path(proj_dir, "figures_manuscript")
source_dir <- file.path(out_dir, "source_data")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)

pal <- c(ink = "#222222", blue = "#2166AC", red = "#B2182B",
         light_blue = "#D9E8F5", grey = "#777777", light_grey = "#E6E6E6")

theme_pub <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = base_size - 0.3, colour = "black"),
      axis.title = element_text(size = base_size),
      strip.text = element_text(size = base_size, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = pal["grey"]),
      plot.caption = element_text(size = base_size - 1, colour = pal["grey"], hjust = 0),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.margin = margin(6, 8, 6, 8)
    )
}
theme_set(theme_pub())

save_pub <- function(plot, stem, width_mm = 183, height_mm = 115, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svglite(file.path(out_dir, paste0(stem, ".svg")), width = w, height = h)
  print(plot); dev.off()
  cairo_pdf(file.path(out_dir, paste0(stem, ".pdf")), width = w, height = h, family = "Arial")
  print(plot); dev.off()
  agg_tiff(file.path(out_dir, paste0(stem, ".tiff")), width = w, height = h,
           units = "in", res = dpi, compression = "lzw")
  print(plot); dev.off()
  agg_png(file.path(out_dir, paste0(stem, ".png")), width = w, height = h,
          units = "in", res = 300)
  print(plot); dev.off()
}

# Figure 1: study flow and estimand definition --------------------------------
flow <- read.csv(file.path(proj_dir, "output", "first_return_selection_flow.csv")) %>%
  mutate(
    study = recode(study, China = "CHARLS", USA = "HRS", England = "ELSA", Europe = "SHARE"),
    retained_label = paste0(format(repeated_observers, big.mark = ","), " repeated observers"),
    eligible_label = paste0(format(baseline_eligible, big.mark = ","), " baseline eligible"),
    excluded_label = paste0(format(singletons, big.mark = ","), " with one assessment"),
    retention_label = paste0(sprintf("%.1f", followup_percent), "% retained"),
    study = factor(study, levels = c("CHARLS", "HRS", "ELSA", "SHARE"))
  )
write.csv(flow, file.path(source_dir, "Figure1_study_flow.csv"), row.names = FALSE)

flow_long <- flow %>%
  select(study, baseline_eligible, repeated_observers) %>%
  pivot_longer(-study, names_to = "stage", values_to = "n") %>%
  mutate(
    stage = factor(stage, levels = c("baseline_eligible", "repeated_observers"),
                   labels = c("Baseline eligible", "Repeated observers")),
    study = factor(study, levels = c("CHARLS", "HRS", "ELSA", "SHARE"))
  )
p1a <- ggplot(flow_long, aes(stage, n, group = study, colour = study)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  geom_text(aes(label = format(n, big.mark = ",")), size = 2.2,
            vjust = -0.8, show.legend = FALSE) +
  geom_text(data = flow %>% mutate(stage = "Repeated observers", n = repeated_observers),
            aes(stage, n, label = retention_label), colour = "grey35",
            size = 2.0, vjust = 1.8, inherit.aes = FALSE) +
  facet_wrap(~study, scales = "free_y", nrow = 1) +
  scale_colour_manual(values = c(CHARLS = unname(pal["red"]), HRS = unname(pal["blue"]),
                                 ELSA = "#4D9221", SHARE = "#762A83")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.22)), labels = scales::comma) +
  labs(x = NULL, y = "Participants") +
  theme_pub(6.6) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        strip.background = element_blank())

estimand_df <- tibble(
  x = c(0, 1, 2, 3),
  label = c("Baseline\nassessment", "First-return\nselection weight",
            "Conditional-response\nIPCW", "Mixed-effects\ntrajectory model")
)
p1b <- ggplot(estimand_df, aes(x, 0)) +
  annotate("segment", x = 0, xend = 3, y = 0, yend = 0, linewidth = 0.7,
           colour = pal["grey"], arrow = arrow(length = unit(2, "mm"))) +
  geom_point(size = 5.2, shape = 21, stroke = 0.6, fill = pal["light_blue"], colour = pal["blue"]) +
  geom_text(aes(label = label), y = -0.22, size = 2.3, lineheight = 0.92) +
  annotate("label", x = 1.5, y = 0.42,
           label = "Primary: Mid−Low difference in average annual memory change, 0–4 years",
           size = 2.5, linewidth = 0.25, fill = "white") +
  annotate("text", x = 1.5, y = 0.23,
           label = "Secondary: difference in instantaneous rate at year 4",
           size = 2.25, colour = pal["grey"]) +
  coord_cartesian(xlim = c(-0.25, 3.25), ylim = c(-0.48, 0.58), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(plot.margin = margin(4, 12, 8, 12))

fig1 <- p1a / p1b + plot_layout(heights = c(1.1, 0.9)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 8, face = "bold"))
save_pub(fig1, "Figure1_study_flow_and_estimand", height_mm = 118)

# Figure 2: primary within-cohort estimates -----------------------------------
primary <- read.csv(file.path(proj_dir, "output", "final_revised_primary_results.csv"))
within <- primary %>%
  filter(estimand == "average_annual_change_difference_0_to_4",
         comparison_type == "within_cohort") %>%
  mutate(
    cohort = recode(comparison, China = "CHARLS", USA = "HRS",
                    England = "ELSA", Europe = "Pooled SHARE"),
    cohort = factor(cohort, levels = rev(c("CHARLS", "HRS", "ELSA", "Pooled SHARE")))
  )
write.csv(within, file.path(source_dir, "Figure2_primary_within_cohort.csv"), row.names = FALSE)
p2 <- ggplot(within, aes(estimate, cohort)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = pal["grey"]) +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), orientation = "y",
                width = 0.18, linewidth = 0.55, colour = pal["ink"]) +
  geom_point(aes(fill = cohort == "CHARLS"), shape = 21, size = 3,
             stroke = 0.55, colour = pal["ink"]) +
  scale_fill_manual(values = c(`TRUE` = unname(pal["red"]), `FALSE` = "white")) +
  scale_x_continuous(breaks = seq(-0.05, 0.15, 0.05)) +
  labs(
    x = "Intermediate minus low education difference\nin average annual memory change, 0–4 years (points/year)",
    y = NULL,
    caption = "Positive values favour intermediate education; error bars show 95% CIs."
  ) + theme_pub(7.3)
save_pub(p2, "Figure2_primary_within_cohort", width_mm = 120, height_mm = 88)

# Figure 3: formal cross-cohort differences -----------------------------------
clean_comp <- function(x) {
  case_when(
    grepl("USA|HRS", x) ~ "CHARLS minus HRS",
    grepl("England|ELSA", x) ~ "CHARLS minus ELSA",
    grepl("Europe|SHARE", x) ~ "CHARLS minus pooled SHARE",
    TRUE ~ x
  )
}
cross <- primary %>%
  filter(comparison_type == "cross_cohort") %>%
  mutate(
    panel = if_else(grepl("average", estimand), "Average annual change, 0–4 years",
                    "Instantaneous rate at year 4"),
    label = clean_comp(comparison),
    label = factor(label, levels = rev(c("CHARLS minus HRS", "CHARLS minus ELSA",
                                        "CHARLS minus pooled SHARE")))
  )
write.csv(cross, file.path(source_dir, "Figure3_cross_cohort_comparisons.csv"), row.names = FALSE)
p3 <- ggplot(cross, aes(estimate, label)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = pal["grey"]) +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), orientation = "y",
                width = 0.18, linewidth = 0.55, colour = pal["ink"]) +
  geom_point(size = 2.8, shape = 21, stroke = 0.55, fill = pal["red"], colour = pal["ink"]) +
  facet_wrap(~panel, nrow = 1, scales = "free_x") +
  labs(
    x = "Difference between education contrasts (points/year)", y = NULL,
    caption = "Positive estimates indicate a larger intermediate-versus-low contrast in CHARLS. Error bars are 95% CIs."
  ) + theme_pub(7) +
  theme(strip.background = element_blank(), strip.text = element_text(hjust = 0))
save_pub(p3, "Figure3_cross_cohort_comparisons", height_mm = 88)

# Figure 4: SHARE country heterogeneity ---------------------------------------
country <- read.csv(file.path(proj_dir, "output", "share_country_common_model_results.csv")) %>%
  filter(estimand == "average_annual_change_difference_0_to_4",
         status %in% c("ok", "singular_random_slope")) %>%
  arrange(estimate) %>%
  mutate(country = factor(country, levels = country))
meta <- read.csv(file.path(proj_dir, "output", "share_country_random_effects_summary.csv")) %>%
  filter(estimand == "average_annual_change_difference_0_to_4")
charls <- within %>% filter(cohort == "CHARLS") %>% pull(estimate)
write.csv(country, file.path(source_dir, "Figure4_SHARE_country_estimates.csv"), row.names = FALSE)
write.csv(meta, file.path(source_dir, "Figure4_SHARE_random_effects_summary.csv"), row.names = FALSE)
p4 <- ggplot(country, aes(estimate, country)) +
  annotate("rect", xmin = meta$prediction_lo, xmax = meta$prediction_hi,
           ymin = -Inf, ymax = Inf, fill = pal["light_blue"], alpha = 0.45) +
  geom_vline(xintercept = 0, colour = pal["grey"], linewidth = 0.4) +
  geom_vline(xintercept = charls, colour = pal["red"], linetype = "dashed", linewidth = 0.65) +
  geom_errorbar(aes(xmin = CI_lo, xmax = CI_hi), orientation = "y",
                width = 0.17, linewidth = 0.45, colour = pal["grey"]) +
  geom_point(aes(shape = singular), size = 2.1, stroke = 0.5,
             fill = "white", colour = pal["blue"]) +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 24)) +
  coord_cartesian(xlim = c(-0.23, 0.14)) +
  labs(
    x = "Intermediate minus low education difference\nin average annual memory change, 0–4 years (points/year)",
    y = NULL,
    caption = paste0(
      "Blue band: 95% random-effects prediction interval; dashed red line: CHARLS estimate (",
      sprintf("%.3f", charls), "). Triangle: boundary-singular random-slope fit."
    )
  ) + theme_pub(6.8)
save_pub(p4, "Figure4_SHARE_country_heterogeneity", width_mm = 150, height_mm = 150)

cat("Manuscript figures written to", out_dir, "\n")
