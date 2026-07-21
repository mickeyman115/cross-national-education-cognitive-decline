###############################################################################
# Alzheimer’s & Dementia submission figures
# Uses locked models/results only. No model fitting is performed.
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(lme4)
  library(svglite)
  library(ragg)
  library(scales)
})

source(if (file.exists("analysis/00_config.R")) "analysis/00_config.R" else "00_config.R")
project_dir <- PROJECT_DIR
figure_dir <- file.path(project_dir, "figures", "ad_submission")
source_dir <- file.path(figure_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

pal <- c(
  ink = "#263238",
  neutral = "#7B8794",
  pale = "#DCE3E8",
  blue = "#2878B5",
  blue_light = "#8EC0DC",
  teal = "#2A9D8F",
  orange = "#E28E2C",
  red = "#C44E52",
  paper = "#FFFFFF"
)

study_cols <- c(
  "CHARLS" = unname(pal["red"]),
  "HRS" = unname(pal["blue"]),
  "ELSA" = unname(pal["teal"]),
  "Pooled SHARE" = unname(pal["orange"])
)

edu_cols <- c(
  "Low" = unname(pal["neutral"]),
  "Intermediate" = unname(pal["blue"]),
  "High" = unname(pal["teal"])
)

model_cols <- c("M0" = unname(pal["neutral"]), "M1" = unname(pal["blue"]))

theme_ad <- function(base_size = 7.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = pal["ink"]),
      axis.ticks = element_line(linewidth = 0.35, colour = pal["ink"]),
      axis.text = element_text(colour = pal["ink"], size = base_size - 0.3),
      axis.title = element_text(colour = pal["ink"], size = base_size),
      legend.title = element_text(size = base_size - 0.2, face = "bold"),
      legend.text = element_text(size = base_size - 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size + 0.2, face = "bold", colour = pal["ink"]),
      plot.title = element_text(size = base_size + 0.8, face = "bold", colour = pal["ink"], hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.1, colour = pal["neutral"], lineheight = 1.05),
      plot.tag = element_text(size = 9, face = "bold", colour = pal["ink"]),
      plot.tag.position = c(0, 1),
      plot.margin = margin(5, 7, 5, 7),
      panel.grid = element_blank(),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(6, "mm")
    )
}

theme_set(theme_ad())

save_publication_figure <- function(plot, stem, width_mm = 183, height_mm = 110, dpi = 600) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  svg_file <- file.path(figure_dir, paste0(stem, ".svg"))
  pdf_file <- file.path(figure_dir, paste0(stem, ".pdf"))
  tif_file <- file.path(figure_dir, paste0(stem, ".tiff"))
  png_file <- file.path(figure_dir, paste0(stem, ".png"))

  svglite::svglite(svg_file, width = w, height = h)
  print(plot)
  dev.off()

  grDevices::cairo_pdf(pdf_file, width = w, height = h, family = "Arial")
  print(plot)
  dev.off()

  ragg::agg_tiff(tif_file, width = w, height = h, units = "in", res = dpi,
                 compression = "lzw", background = "white")
  print(plot)
  dev.off()

  ragg::agg_png(png_file, width = w, height = h, units = "in", res = 300,
                background = "white")
  print(plot)
  dev.off()
}

fmt_n <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
rename_study <- function(x) {
  recode(as.character(x), China = "CHARLS", USA = "HRS", England = "ELSA",
         Europe = "Pooled SHARE", .default = as.character(x))
}

###############################################################################
# Figure 1: design and estimand
###############################################################################

table1 <- read.csv(file.path(project_dir, "docs", "AD_TABLE1_DRAFT.csv"), check.names = FALSE)
cohort_cards <- data.frame(
  cohort = c("CHARLS", "HRS", "ELSA", "Pooled SHARE"),
  context = c("China", "United States", "England", "SHARE countries"),
  participants = as.numeric(table1[table1$Characteristic == "Participants", -1]),
  observations = as.numeric(table1[table1$Characteristic == "Memory observations", -1]),
  x = c(1, 2, 1, 2), y = c(2, 2, 1, 1)
)

p1a <- ggplot(cohort_cards) +
  geom_tile(aes(x, y, fill = cohort), width = 0.88, height = 0.78, colour = "white", linewidth = 0.6) +
  geom_text(aes(x, y + 0.18, label = cohort), colour = "white", family = "Arial", fontface = "bold", size = 3.2) +
  geom_text(aes(x, y - 0.03, label = context), colour = "white", family = "Arial", size = 2.35) +
  geom_text(aes(x, y - 0.23, label = paste0(fmt_n(participants), " participants\n", fmt_n(observations), " observations")),
            colour = "white", family = "Arial", size = 2.05, lineheight = 0.95) +
  scale_fill_manual(values = study_cols) +
  coord_cartesian(xlim = c(0.45, 2.55), ylim = c(0.48, 2.52), clip = "off") +
  labs(title = "Four harmonized aging cohorts") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold", colour = pal["ink"], hjust = 0),
        legend.position = "none", plot.margin = margin(8, 5, 4, 5))

t <- seq(0, 7, length.out = 141)
schematic <- bind_rows(
  data.frame(time = t, score = 10.2 - 0.10 * t - 0.030 * t^2, education = "Low"),
  data.frame(time = t, score = 10.2 - 0.08 * t - 0.014 * t^2, education = "Intermediate")
)
tangent <- data.frame(
  education = c("Low", "Intermediate"),
  slope = c(-0.10 - 2 * 0.030 * 4, -0.08 - 2 * 0.014 * 4),
  y4 = c(10.2 - 0.10 * 4 - 0.030 * 16, 10.2 - 0.08 * 4 - 0.014 * 16)
) %>%
  rowwise() %>%
  do(data.frame(education = .$education, time = c(3.25, 4.75),
                score = .$y4 + .$slope * (c(3.25, 4.75) - 4))) %>%
  ungroup()

p1b <- ggplot(schematic, aes(time, score, colour = education)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = 4, linetype = "dotted", linewidth = 0.45, colour = pal["ink"]) +
  geom_line(data = tangent, linewidth = 0.7, linetype = "longdash") +
  annotate("text", x = 4.08, y = 10.35, label = "Year 4", hjust = 0, size = 2.4, family = "Arial") +
  annotate("label", x = 4.2, y = 8.55,
           label = "Year-4 rate contrast\n(intermediate − low)",
           size = 2.3, family = "Arial", linewidth = 0.2, fill = "white") +
  scale_colour_manual(values = edu_cols[c("Low", "Intermediate")]) +
  scale_x_continuous(breaks = c(0, 2, 4, 6), limits = c(0, 7)) +
  labs(title = "A rate contrast, not four-year change", x = "Years since cohort entry", y = "Episodic memory score", colour = "Education") +
  theme_ad(7.2) +
  theme(legend.position = "bottom")

flow <- data.frame(
  x = c(1, 2.5, 4.0, 5.5), y = c(1.55, 1.55, 1.55, 1.55),
  label = c("M0\nBase + IPCW", "M1\n+ childhood SEP", "M2 / M3\n+ wealth or income", "M4\n+ wealth + income"),
  type = c("Primary", "Primary", "Exploratory", "Exploratory")
)
flow_edges <- data.frame(
  x = flow$x[-nrow(flow)], xend = flow$x[-1],
  y = flow$y[-nrow(flow)], yend = flow$y[-1]
)
p1c <- ggplot(flow) +
  geom_segment(aes(x = x, xend = xend, y = y, yend = yend),
               data = flow_edges, arrow = arrow(length = unit(1.8, "mm")),
               linewidth = 0.45, colour = pal["neutral"]) +
  geom_label(aes(x, y, label = label, fill = type), size = 2.35, family = "Arial",
             fontface = "bold", linewidth = 0.25, label.padding = unit(1.8, "mm"), colour = pal["ink"]) +
  annotate("text", x = 1.75, y = 0.72,
           label = "Full sample: 144,642 participants / 542,426 observations",
           family = "Arial", size = 2.3, colour = pal["ink"]) +
  annotate("text", x = 4.75, y = 0.72,
           label = "Common economic sample: 131,588 / 511,722",
           family = "Arial", size = 2.3, colour = pal["ink"]) +
  scale_fill_manual(values = c(Primary = unname(pal["blue_light"]), Exploratory = "#F5D4A6")) +
  coord_cartesian(xlim = c(0.35, 6.15), ylim = c(0.35, 2.15), clip = "off") +
  labs(title = "Analysis hierarchy") +
  theme_void(base_family = "Arial") +
  theme(plot.title = element_text(size = 8, face = "bold", colour = pal["ink"]), legend.position = "none")

fig1 <- (p1a | p1b) / p1c +
  plot_layout(heights = c(1.65, 1), widths = c(1, 1.25)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))

save_publication_figure(fig1, "Figure1_study_design_estimand", 183, 116)
write.csv(cohort_cards[, c("cohort", "context", "participants", "observations")],
          file.path(source_dir, "Figure1_cohort_counts.csv"), row.names = FALSE)

###############################################################################
# Figure 2: M0 quadratic trajectories from locked model
###############################################################################

m0 <- readRDS(file.path(project_dir, "output", "fullsample_m0_model.rds"))
mf <- model.frame(m0)
support <- mf %>%
  transmute(study = as.character(study), education = as.character(edu3_f),
            t_round = round(time_in_study)) %>%
  count(study, education, t_round, name = "n_obs")

max_support <- support %>%
  filter(n_obs >= 50) %>%
  group_by(study, education) %>%
  summarise(max_time = max(t_round), .groups = "drop")

prediction_grid <- max_support %>%
  rowwise() %>%
  do(data.frame(study = .$study, education = .$education,
                time_in_study = seq(0, .$max_time, by = 0.1))) %>%
  ungroup() %>%
  mutate(
    age_base_c = 0,
    enroll_year_c = 0,
    retest_flag = 1,
    female = 0.5,
    study = factor(study, levels = levels(mf$study)),
    edu3_f = factor(education, levels = levels(mf$edu3_f))
  )

fixed_formula <- lme4::nobars(formula(m0))
X <- model.matrix(delete.response(terms(fixed_formula)), prediction_grid)
beta <- lme4::fixef(m0)
X <- X[, names(beta), drop = FALSE]
V <- as.matrix(vcov(m0))
prediction_grid$estimate <- as.numeric(X %*% beta)
prediction_grid$SE <- sqrt(rowSums((X %*% V) * X))
prediction_grid$CI_lo <- prediction_grid$estimate - 1.96 * prediction_grid$SE
prediction_grid$CI_hi <- prediction_grid$estimate + 1.96 * prediction_grid$SE
prediction_grid$study_label <- factor(rename_study(prediction_grid$study),
                                      levels = c("CHARLS", "HRS", "ELSA", "Pooled SHARE"))
prediction_grid$education <- factor(recode(prediction_grid$education, Mid = "Intermediate"),
                                    levels = c("Low", "Intermediate", "High"))

write.csv(prediction_grid %>%
            transmute(cohort = study_label, education, years_since_entry = time_in_study,
                      predicted_score = estimate, SE, CI_lo, CI_hi),
          file.path(source_dir, "Figure2_model_trajectories.csv"), row.names = FALSE)
write.csv(max_support %>% mutate(study = rename_study(study), education = recode(education, Mid = "Intermediate")),
          file.path(source_dir, "Figure2_support_limits.csv"), row.names = FALSE)

p2 <- ggplot(prediction_grid,
             aes(time_in_study, estimate, colour = education, fill = education,
                 linetype = education, group = education)) +
  geom_ribbon(aes(ymin = CI_lo, ymax = CI_hi), alpha = 0.11, colour = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 4, linewidth = 0.38, linetype = "dotted", colour = pal["ink"]) +
  facet_wrap(~ study_label, nrow = 1, scales = "free_x") +
  scale_colour_manual(values = edu_cols) +
  scale_fill_manual(values = edu_cols) +
  scale_linetype_manual(values = c(Low = "solid", Intermediate = "solid", High = "22")) +
  scale_x_continuous(breaks = seq(0, 16, 2), expand = expansion(mult = c(0.01, 0.02))) +
  labs(x = "Years since cohort entry", y = "Predicted episodic memory score (0–20)",
       colour = "Education", fill = "Education", linetype = "Education") +
  theme_ad(7.2) +
  theme(legend.position = "bottom", legend.box = "horizontal",
        panel.spacing.x = unit(4, "mm"))

save_publication_figure(p2, "Figure2_quadratic_memory_trajectories", 183, 90)

###############################################################################
# Figure 3: primary forest plots
###############################################################################

primary <- read.csv(file.path(project_dir, "output", "final_primary_results.csv"), check.names = FALSE) %>%
  mutate(model = factor(model, levels = c("M0", "M1")))

within <- primary %>%
  filter(contrast == "Within Mid-Low") %>%
  mutate(label = rename_study(study),
         label = factor(label, levels = rev(c("CHARLS", "HRS", "ELSA", "Pooled SHARE"))),
         y_base = as.numeric(label),
         y = y_base + ifelse(model == "M0", -0.12, 0.12),
         value = sprintf("%.3f (%.3f, %.3f)", estimate, CI_lo, CI_hi))

cross <- primary %>%
  filter(contrast == "Cross-national") %>%
  mutate(label = recode(study,
                        "China vs USA" = "CHARLS − HRS",
                        "China vs England" = "CHARLS − ELSA",
                        "China vs Europe" = "CHARLS − pooled SHARE"),
         label = factor(label, levels = rev(c("CHARLS − HRS", "CHARLS − ELSA", "CHARLS − pooled SHARE"))),
         y_base = as.numeric(label),
         y = y_base + ifelse(model == "M0", -0.12, 0.12),
         value = sprintf("%.3f (%.3f, %.3f)", estimate, CI_lo, CI_hi))

write.csv(within %>% select(model, cohort = label, estimate, SE, CI_lo, CI_hi, p_value, FMI, MC_error),
          file.path(source_dir, "Figure3a_within_cohort_contrasts.csv"), row.names = FALSE)
write.csv(cross %>% select(model, comparison = label, estimate, SE, CI_lo, CI_hi, p_value, p_holm, FMI, MC_error),
          file.path(source_dir, "Figure3b_cross_cohort_contrasts.csv"), row.names = FALSE)

p3a <- ggplot(within, aes(colour = model)) +
  geom_vline(xintercept = 0, linewidth = 0.38, linetype = "dashed", colour = pal["neutral"]) +
  geom_segment(aes(x = CI_lo, xend = CI_hi, y = y, yend = y), linewidth = 0.7) +
  geom_point(aes(x = estimate, y = y, shape = model), size = 2.2, stroke = 0.5) +
  scale_colour_manual(values = model_cols) +
  scale_shape_manual(values = c(M0 = 16, M1 = 15)) +
  scale_y_continuous(breaks = sort(unique(within$y_base)),
                     labels = levels(within$label), expand = expansion(add = c(0.45, 0.45))) +
  scale_x_continuous(limits = c(-0.04, 0.28), breaks = seq(-0.05, 0.30, 0.05)) +
  labs(title = "Within-cohort contrasts", x = "Intermediate − low difference in year-4 rate\n(memory points/year)", y = NULL,
       colour = "Model", shape = "Model") +
  theme_ad(7.2) +
  theme(legend.position = "bottom")

p3b <- ggplot(cross, aes(colour = model)) +
  geom_vline(xintercept = 0, linewidth = 0.38, linetype = "dashed", colour = pal["neutral"]) +
  geom_segment(aes(x = CI_lo, xend = CI_hi, y = y, yend = y), linewidth = 0.7) +
  geom_point(aes(x = estimate, y = y, shape = model), size = 2.2, stroke = 0.5) +
  scale_colour_manual(values = model_cols) +
  scale_shape_manual(values = c(M0 = 16, M1 = 15)) +
  scale_y_continuous(breaks = sort(unique(cross$y_base)),
                     labels = levels(cross$label), expand = expansion(add = c(0.45, 0.45))) +
  scale_x_continuous(limits = c(0, 0.29), breaks = seq(0, 0.30, 0.05)) +
  labs(title = "Formal cross-cohort differences", x = "Difference between education contrasts\n(memory points/year)", y = NULL,
       colour = "Model", shape = "Model") +
  theme_ad(7.2) +
  theme(legend.position = "bottom")

fig3 <- (p3a | p3b) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = "bottom", plot.tag = element_text(size = 9, face = "bold", family = "Arial"))

save_publication_figure(fig3, "Figure3_primary_contrasts", 183, 98)

###############################################################################
# Figure 4: sequential adjustment and sensitivity analyses
###############################################################################

mechanism <- read.csv(file.path(project_dir, "output", "final_mechanism_results.csv"), check.names = FALSE) %>%
  filter(study == "China", contrast == "Mid - Low") %>%
  mutate(
    model_label = recode(model,
                         M0_Base = "M0  Base",
                         M1_ChildSES = "M1  + childhood SEP",
                         M2_Wealth = "M2  + wealth",
                         M3_Income = "M3  + income",
                         M4_Joint = "M4  + wealth + income"),
    model_label = factor(model_label, levels = rev(c("M0  Base", "M1  + childhood SEP", "M2  + wealth", "M3  + income", "M4  + wealth + income"))),
    attenuation = 100 * (Q_bar[model == "M0_Base"] - Q_bar) / Q_bar[model == "M0_Base"],
    attenuation_label = ifelse(model == "M0_Base", "Reference", sprintf("−%.1f%%", attenuation)),
    family = case_when(model == "M0_Base" ~ "Base", model == "M1_ChildSES" ~ "Childhood SEP", TRUE ~ "Later-life economics")
  )

write.csv(mechanism %>% select(model, estimate = Q_bar, SE = SE_pooled, CI_lo, CI_hi, attenuation, FMI, MC_error),
          file.path(source_dir, "Figure4a_sequential_adjustment.csv"), row.names = FALSE)

p4a <- ggplot(mechanism, aes(Q_bar, model_label, colour = family)) +
  geom_vline(xintercept = mechanism$Q_bar[mechanism$model == "M0_Base"], linetype = "dotted",
             linewidth = 0.4, colour = pal["neutral"]) +
  geom_segment(aes(x = CI_lo, xend = CI_hi, y = model_label, yend = model_label), linewidth = 0.75) +
  geom_point(size = 2.4) +
  geom_text(aes(x = 0.274, label = attenuation_label), hjust = 1, size = 2.3, colour = pal["ink"], family = "Arial") +
  scale_colour_manual(values = c("Base" = unname(pal["neutral"]),
                                 "Childhood SEP" = unname(pal["blue"]),
                                 "Later-life economics" = unname(pal["orange"]))) +
  scale_x_continuous(limits = c(0.15, 0.278), breaks = seq(0.16, 0.28, 0.02)) +
  labs(title = "Sequential socioeconomic adjustment", subtitle = "Common economic-data sample (n = 131,588)",
       x = "CHARLS intermediate − low contrast\n(memory points/year)", y = NULL, colour = NULL) +
  theme_ad(7.2) +
  theme(legend.position = "bottom")

ipcw <- read.csv(file.path(project_dir, "output", "ipcw_sensitivity.csv"), check.names = FALSE) %>%
  filter(study == "China", contrast == "Within Mid-Low") %>%
  mutate(specification = recode(model,
                                IPCW_current = "Primary IPCW",
                                Unweighted = "No attrition weight",
                                IPCW_truncated = "Alternate truncation"),
         specification = factor(specification, levels = rev(c("Primary IPCW", "No attrition weight", "Alternate truncation"))))

write.csv(ipcw %>% select(specification, estimate, SE, CI_lo, CI_hi),
          file.path(source_dir, "Figure4b_ipcw_sensitivity.csv"), row.names = FALSE)

p4b <- ggplot(ipcw, aes(estimate, specification)) +
  geom_vline(xintercept = ipcw$estimate[ipcw$specification == "Primary IPCW"], linetype = "dotted",
             linewidth = 0.4, colour = pal["neutral"]) +
  geom_segment(aes(x = CI_lo, xend = CI_hi, y = specification, yend = specification), colour = pal["blue"], linewidth = 0.75) +
  geom_point(colour = pal["blue"], size = 2.3) +
  scale_x_continuous(limits = c(0.165, 0.268), breaks = seq(0.17, 0.27, 0.02)) +
  labs(title = "Attrition-weight specification", x = "CHARLS contrast (memory points/year)", y = NULL) +
  theme_ad(7.2)

charls_wave <- read.csv(file.path(project_dir, "output", "charls_sensitivity.csv"), check.names = FALSE) %>%
  filter(wave == "Wave_2011", model %in% c("M0", "M1")) %>%
  transmute(sample = "2011-entry cohort", model, estimate, SE, CI_lo, CI_hi)
charls_full <- primary %>%
  filter(study == "China", contrast == "Within Mid-Low", model %in% c("M0", "M1")) %>%
  transmute(sample = "Full CHARLS sample", model = as.character(model), estimate, SE, CI_lo, CI_hi)
charls_compare <- bind_rows(charls_full, charls_wave) %>%
  mutate(model = factor(model, levels = c("M0", "M1")),
         sample = factor(sample, levels = rev(c("Full CHARLS sample", "2011-entry cohort"))),
         y_base = as.numeric(sample), y = y_base + ifelse(model == "M0", -0.10, 0.10))

write.csv(charls_compare %>% select(sample, model, estimate, SE, CI_lo, CI_hi),
          file.path(source_dir, "Figure4c_charls_entry_sensitivity.csv"), row.names = FALSE)

p4c <- ggplot(charls_compare, aes(colour = model)) +
  geom_segment(aes(x = CI_lo, xend = CI_hi, y = y, yend = y), linewidth = 0.7) +
  geom_point(aes(x = estimate, y = y, shape = model), size = 2.2) +
  scale_colour_manual(values = model_cols) +
  scale_shape_manual(values = c(M0 = 16, M1 = 15)) +
  scale_y_continuous(breaks = sort(unique(charls_compare$y_base)), labels = levels(charls_compare$sample),
                     expand = expansion(add = c(0.4, 0.4))) +
  scale_x_continuous(limits = c(0.165, 0.30), breaks = seq(0.17, 0.29, 0.03)) +
  labs(title = "CHARLS entry-wave sensitivity", x = "Intermediate − low contrast (memory points/year)", y = NULL,
       colour = "Model", shape = "Model") +
  theme_ad(7.2) +
  theme(legend.position = "bottom")

fig4 <- (p4a | (p4b / p4c)) +
  plot_layout(widths = c(1.12, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 9, face = "bold", family = "Arial"))

save_publication_figure(fig4, "Figure4_adjustment_and_sensitivity", 183, 122)

###############################################################################
# QA inventory
###############################################################################

inventory <- data.frame(
  figure = paste0("Figure", 1:4),
  stem = c("Figure1_study_design_estimand", "Figure2_quadratic_memory_trajectories",
           "Figure3_primary_contrasts", "Figure4_adjustment_and_sensitivity"),
  width_mm = 183,
  height_mm = c(116, 90, 98, 122),
  formats = "SVG; PDF; TIFF 600 dpi; PNG 300 dpi",
  stringsAsFactors = FALSE
)
write.csv(inventory, file.path(figure_dir, "figure_inventory.csv"), row.names = FALSE)

cat("Completed four Alzheimer’s & Dementia submission figures in:\n", figure_dir, "\n")
