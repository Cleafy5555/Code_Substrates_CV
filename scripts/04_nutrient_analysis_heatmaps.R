## ================================================================
## Nutrient analysis + heatmaps
## ================================================================

## ---- Libraries ----
library(tidyverse)
library(readxl)
library(janitor)
library(stringr)
library(forcats)
library(rstatix)
library(emmeans)
library(broom)
library(writexl)
library(scales)

## ---- Paths ----
root    <- "//ad.uws.edu.au/dfshare/HomesHWK$/90955975/Desktop/Claudia/PhD OG werk"
gt_root <- file.path(root, "GROWTH TRIAL")

new_file <- "TISSUEa_80503_092625.xls"
old_file <- file.path(gt_root, "broc_nutrient.csv")

out_dir <- file.path(root, "THESIS/Chapter 1/Figures")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

xlsx_out          <- file.path(out_dir, "additional_nutrient_stats_table.xlsx")
png_out_facet     <- file.path(out_dir, "additional_nutrient_heatmap_facets.png")
svg_out_facet     <- file.path(out_dir, "additional_nutrient_heatmap_facets.svg")
png_out_onepanel  <- file.path(out_dir, "additional_nutrient_heatmap_onepanel_PUBLICATION.png")
svg_out_onepanel  <- file.path(out_dir, "additional_nutrient_heatmap_onepanel_PUBLICATION.svg")

## ================================================================
## Helpers
## ================================================================

## Robust numeric parsing (handles "<0.5", "0,45", ">10", etc.)
parse_num <- function(x) readr::parse_number(as.character(x))

## Common nutrient column detector
nutrient_pattern <- paste0(
  "(",
  paste(c(
    "^n(_percent)?$", "^p(_percent)?$", "^k(_percent)?$", "^ca(_percent)?$",
    "^mg(_percent)?$", "^s(_percent)?$",
    "boron|^b(_|$)", "iron|^fe(_|$)", "zinc|^zn(_|$)", "manganese|^mn(_|$)",
    "copper|^cu(_|$)", "sodium|^na(_|$)", "potassium|^k(_|$)", "calcium|^ca(_|$)",
    "magnesium|^mg(_|$)", "sulfur|sulphur|^s(_|$)",
    "alumini?um|^al(_|$)"
  ), collapse = "|"),
  ")"
)

## Pairwise Wilcoxon vs. control (BH) → p + significance code
pairwise_vs_control <- function(df) {
  ctrl <- levels(df$treatment)[1]
  res <- rstatix::wilcox_test(
    df, value ~ treatment,
    ref.group = ctrl, p.adjust.method = "BH"
  ) %>%
    rstatix::add_significance(p.col = "p.adj")
  
  if ("p.adj.signif" %in% names(res)) {
    res <- dplyr::rename(res, sig = p.adj.signif)
  } else if ("p.signif" %in% names(res)) {
    res <- dplyr::rename(res, sig = p.signif)
  } else {
    res <- dplyr::mutate(res, sig = NA_character_)
  }
  
  res %>%
    mutate(treatment = group2) %>%
    select(treatment, p = p.adj, sig)
}

## Core analysis for any run
analyze_run <- function(df_wide,
                        run_label       = "Run_new",
                        tissue_filter   = c("stem", "stems"),
                        control_label   = "control",
                        treatment_levels = NULL) {
  # Ensure minimal columns exist
  if (!"tissue"   %in% names(df_wide)) df_wide$tissue   <- factor("stem")
  if (!"infected" %in% names(df_wide)) df_wide$infected <- factor("unknown")
  if (!"id"       %in% names(df_wide)) df_wide$id       <- seq_len(nrow(df_wide))
  
  # Filter tissue
  df_wide <- df_wide %>%
    mutate(tissue = as.character(tissue)) %>%
    filter(str_to_lower(tissue) %in% str_to_lower(tissue_filter)) %>%
    mutate(tissue = factor(tissue))
  
  # Parse numeric nutrient columns
  df_wide <- df_wide %>%
    mutate(across(matches(nutrient_pattern), parse_num, .names = "{.col}"))
  
  meta_cols <- c("id", "run", "treatment", "infected", "tissue")
  num_cols  <- names(df_wide)[sapply(df_wide, is.numeric)]
  nutrient_cols <- setdiff(num_cols, intersect(num_cols, meta_cols))
  stopifnot(length(nutrient_cols) > 0)
  
  # Treatment factoring
  if (!is.null(treatment_levels)) {
    df_wide <- df_wide %>%
      mutate(treatment = factor(as.character(treatment),
                                levels = treatment_levels)) %>%
      droplevels()
  }
  # Force control first
  df_wide <- df_wide %>%
    mutate(treatment = fct_relevel(factor(treatment), control_label))
  
  # Long format
  long <- df_wide %>%
    pivot_longer(all_of(nutrient_cols), names_to = "nutrient", values_to = "value") %>%
    filter(!is.na(value)) %>%
    mutate(run = run_label)
  
  # Medians + IQR
  medians <- long %>%
    group_by(nutrient, treatment) %>%
    summarise(
      n      = sum(!is.na(value)),
      median = median(value, na.rm = TRUE),
      iqr    = IQR(value, na.rm = TRUE),
      .groups = "drop"
    )
  
  ctrl_meds <- medians %>%
    filter(treatment == control_label) %>%
    select(nutrient, ctrl_med = median)
  
  summary_tbl <- medians %>%
    left_join(ctrl_meds, by = "nutrient") %>%
    mutate(pct_change_vs_ctrl = 100 * (median - ctrl_med) / ctrl_med)
  
  # Wilcoxon vs control
  wilcoxon_results <- long %>%
    mutate(treatment = fct_relevel(treatment, control_label)) %>%
    group_by(nutrient) %>%
    group_modify(~ pairwise_vs_control(.x)) %>%
    ungroup()
  
  summary_tbl <- summary_tbl %>%
    left_join(wilcoxon_results, by = c("nutrient", "treatment")) %>%
    mutate(
      sig = replace_na(sig, "ns"),
      run = run_label
    ) %>%
    arrange(nutrient, treatment)
  
  list(long = long, summary = summary_tbl)
}

## ================================================================
## NEW dataset (TISSUEa_80503_092625.xls)
## ================================================================
raw_new <- read_excel(file.path(gt_root, new_file)) %>% clean_names()
stopifnot("sample_description_number_2" %in% names(raw_new))

dat_new <- raw_new %>%
  rename(id = sample_description_number_2) %>%
  mutate(
    id    = as.character(id) |> str_trim(),
    id_lc = str_to_upper(id),
    # precedence: Control (C*) > New rockwool (NR*) > Used (suffix S/B)
    treatment = case_when(
      str_starts(id_lc, "C")  ~ "control",
      str_starts(id_lc, "NR") ~ "n_rockwool",
      str_ends(id_lc, "S")    ~ "u_rockwools",  # tomato
      str_ends(id_lc, "B")    ~ "u_rockwoolb",  # cucumber
      TRUE                    ~ "unknown"
    ),
    treatment = factor(
      treatment,
      levels = c("control", "n_rockwool", "u_rockwools", "u_rockwoolb", "unknown")
    ),
    run      = "Run_new",
    tissue   = if ("tissue"   %in% names(raw_new)) as.factor(.data$tissue)   else factor("stem"),
    infected = if ("infected" %in% names(raw_new)) as.factor(.data$infected) else factor("unknown")
  ) %>%
  select(-id_lc)

new_res <- analyze_run(
  df_wide         = dat_new,
  run_label       = "Run_new",
  tissue_filter   = c("stem", "stems"),
  control_label   = "control",
  treatment_levels = c("control","n_rockwool","u_rockwools","u_rockwoolb","unknown")
)

## ================================================================
## OLD dataset (broc_nutrient.csv)
## ================================================================
broc_raw <- read.csv(old_file, stringsAsFactors = FALSE)

# Drop all-NA columns
broc_clean <- broc_raw[, colSums(is.na(broc_raw)) < nrow(broc_raw)]

# Map numeric codes → human names if present in 3rd column
if (ncol(broc_clean) >= 3) {
  substrate_names <- c(
    "0"  = "Control",
    "1"  = "Recycled textile",
    "2"  = "Paper waste",
    "3"  = "Perlite",
    "4"  = "Biochar",
    "5"  = "Enriched Biochar",
    "6"  = "Coco coir",
    "7"  = "Sawdust",
    "8"  = "Peat",
    "9"  = "Mushroom waste",
    "10" = "Coffee grounds",
    "11" = "Sand vermicompost",
    "12" = "Peat2"
  )
  
  broc_clean[[3]] <- as.character(broc_clean[[3]])
  code_matches <- broc_clean[[3]] %in% names(substrate_names)
  
  if (any(code_matches)) {
    broc_clean$treatment_name <- broc_clean[[3]]
    broc_clean$treatment_name[code_matches] <-
      substrate_names[broc_clean[[3]][code_matches]]
  } else if ("treatment" %in% tolower(names(broc_clean))) {
    nm_low <- tolower(names(broc_clean))
    broc_clean$treatment_name <- broc_clean[[which(nm_low == "treatment")]]
  } else {
    broc_clean$treatment_name <- broc_clean[[3]]
  }
} else {
  stop("Old dataset has unexpected shape; need at least 3 columns.")
}

# Normalise treatment names
normalize_trt <- function(x) {
  x %>%
    tolower() %>%
    str_replace_all(" ", "_") %>%
    str_replace_all("-", "_")
}

dat_old <- broc_clean %>%
  mutate(
    treatment = normalize_trt(treatment_name),
    tissue    = if ("tissue"   %in% names(.)) .data$tissue   else "stem",
    infected  = if ("infected" %in% names(.)) .data$infected else "unknown",
    run       = "Run_old"
  ) %>%
  { if (!"sample_id" %in% names(.)) mutate(., sample_id = row_number()) else . } %>%
  rename(id = sample_id)

# Coerce any control variants → "control"
dat_old <- dat_old %>%
  mutate(treatment = case_when(
    treatment %in% c("control","0","s0","control_0","_control") ~ "control",
    TRUE ~ treatment
  ))

old_levels <- c(
  "control",
  "recycled_textile","paper_waste","perlite","biochar","enriched_biochar",
  "coco_coir","sawdust","peat","mushroom_waste","coffee_grounds",
  "sand_vermicompost","peat2",
  "n_rockwool","u_rockwools","u_rockwoolb"
)

old_res <- analyze_run(
  df_wide         = dat_old,
  run_label       = "Run_old",
  tissue_filter   = c("stem", "stems"),
  control_label   = "control",
  treatment_levels = old_levels
)

## ================================================================
## Combine & export tables
## ================================================================
summary_all <- bind_rows(new_res$summary, old_res$summary) %>%
  arrange(run, nutrient, treatment)

summary_wide <- summary_all %>%
  mutate(label = sprintf("%.2f (Δ%.0f%%)%s",
                         median, pct_change_vs_ctrl,
                         if_else(sig == "ns", "", paste0(" ", sig)))) %>%
  select(run, nutrient, treatment, label) %>%
  pivot_wider(names_from = c(run, treatment), values_from = label) %>%
  arrange(nutrient)

write_xlsx(
  list(
    "Summary_all_long" = summary_all,
    "Summary_all_wide" = summary_wide,
    "NEW_only_long"    = new_res$summary,
    "OLD_only_long"    = old_res$summary
  ),
  path = xlsx_out
)

message("Saved tables to: ", xlsx_out)

## ================================================================
## Heatmap prep (shared for all plots)
## ================================================================

# base data for heatmaps
hm_df <- summary_all %>%
  select(run, nutrient, treatment, pct_change_vs_ctrl, sig)

# Harmonise nutrient names
nutrient_map <- c(
  "zinc_ppm"           = "Zinc", "zn" = "Zinc", "zn_mg_kg" = "Zinc", "Zinc" = "Zinc",
  "aluminum_ppm"       = "Aluminum", "aluminium_ppm" = "Aluminum",
  "al_mg_kg"           = "Aluminum", "al" = "Aluminum", "Aluminum" = "Aluminum",
  "iron_ppm"           = "Iron", "fe" = "Iron", "Fe" = "Iron", "Iron" = "Iron",
  "copper_ppm"         = "Copper", "cu" = "Copper", "Copper" = "Copper",
  "manganese_ppm"      = "Manganese", "mn" = "Manganese", "Manganese" = "Manganese",
  "boron_ppm"          = "Boron", "b" = "Boron", "Boron" = "Boron",
  "sulfur_percent"     = "Sulfur", "sulphur" = "Sulfur", "Sulfur" = "Sulfur",
  "nitrogen_percent"   = "Nitrogen", "Nitrogen" = "Nitrogen",
  "phosphorous_percent"= "Phosphorus", "Phosphorous" = "Phosphorus",
  "potassium_percent"  = "Potassium", "Potassium" = "Potassium",
  "calcium_percent"    = "Calcium", "Calcium" = "Calcium",
  "magnesium_percent"  = "Magnesium", "Magnesium" = "Magnesium",
  "sodium_ppm"         = "Sodium", "Sodium" = "Sodium"
)

# Plotting orders
treat_order <- c(
  "n_rockwool","u_rockwools","u_rockwoolb",
  "biochar","enriched_biochar","mushroom_waste","sand_vermicompost",
  "recycled_textile","paper_waste","perlite",
  "coco_coir","sawdust","peat","peat2","coffee_grounds","unknown"
)

nutrient_order <- c(
  "Nitrogen","Phosphorus","Potassium","Calcium","Magnesium","Sulfur",
  "Iron","Manganese","Zinc","Copper","Boron","Sodium","Aluminum"
)

# Clean treatment labels for x-axis and legend
legend_labels <- c(
  "control"           = "Control",
  "biochar"           = "Biochar",
  "enriched_biochar"  = "Enriched biochar",
  "mushroom_waste"    = "Mushroom waste",
  "sand_vermicompost" = "Sand vermicompost",
  "u_rockwools"       = "Used rockwool (tomato)",
  "u_rockwoolb"       = "Used rockwool (cucumber)",
  "n_rockwool"        = "New rockwool",
  "recycled_textile"  = "Recycled textile",
  "paper_waste"       = "Paper waste",
  "perlite"           = "Perlite",
  "coco_coir"         = "Coco coir",
  "sawdust"           = "Sawdust",
  "peat"              = "Peat–perlite",
  "peat2"             = "Peat–coir",
  "coffee_grounds"    = "Coffee grounds",
  "unknown"           = "Unknown"
)


# summarise per run × treatment × nutrient (exclude control from plotting)
hm_df_proc <- hm_df %>%
  filter(treatment != "control") %>%
  mutate(
    nutrient  = recode(nutrient, !!!nutrient_map),
    nutrient  = factor(nutrient, levels = rev(nutrient_order)),
    treatment = factor(treatment, levels = treat_order)
  ) %>%
  group_by(run, treatment, nutrient) %>%
  summarise(
    pct_change_vs_ctrl = mean(pct_change_vs_ctrl, na.rm = TRUE),
    sig = paste(unique(sig[sig != "ns"]), collapse = " "),
    .groups = "drop"
  ) %>%
  mutate(sig = ifelse(sig == "", "ns", sig)) %>%
  filter(!is.na(nutrient), !is.na(treatment))

# Cap extremes for colour scale
cap <- 200
hm_df_proc <- hm_df_proc %>%
  mutate(pct_capped = pmax(pmin(pct_change_vs_ctrl, cap), -cap))

## ================================================================
## Heatmap 1 – Faceted by run (diagnostic / supplementary)
## ================================================================
p_hm_facet <- ggplot(hm_df_proc, aes(x = treatment, y = nutrient, fill = pct_capped)) +
  geom_tile() +
  geom_text(aes(label = ifelse(sig == "ns", "", sig)), size = 3.2) +
  scale_fill_gradient2(
    name = "Δ vs control",
    low = muted("blue"), mid = "white", high = muted("red"),
    midpoint = 0, limits = c(-cap, cap),
    labels = function(z) paste0(z, "%")
  ) +
  labs(
    x = NULL, y = NULL,
    title = "Nutrient change vs control (median, Wilcoxon BH significance)",
    subtitle = "Each run shown separately"
  ) +
  facet_wrap(~ run, ncol = 1) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid   = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.text   = element_text(face = "bold")
  )

print(p_hm_facet)

ggsave(png_out_facet, p_hm_facet, width = 12, height = 10, dpi = 300)
ggsave(svg_out_facet, p_hm_facet, width = 12, height = 10)

## ================================================================
## Heatmap 2 – One-panel publication-ready heatmap
## ================================================================

# Build combined x factor with run info
present_new <- hm_df_proc %>%
  filter(run == "Run_new") %>%
  pull(treatment) %>%
  as.character() %>%
  unique()

present_old <- hm_df_proc %>%
  filter(run == "Run_old") %>%
  pull(treatment) %>%
  as.character() %>%
  unique()

present_new <- treat_order[treat_order %in% present_new]
present_old <- treat_order[treat_order %in% present_old]

x_levels <- c(
  paste0("Run_new_", present_new),
  paste0("Run_old_", present_old)
)

hm_df_onepanel <- hm_df_proc %>%
  mutate(
    trt_run      = paste0(run, "_", as.character(treatment)),
    trt_run      = factor(trt_run, levels = x_levels),
    treatment_label = recode(as.character(treatment), !!!legend_labels)
  ) %>%
  droplevels()

divider_x <- length(present_new)  # dashed line between runs

# Colour-blind-friendly diverging palette (print-friendly)
col_low  <- "#3B4CC0"  # deep blue
col_mid  <- "#F7F7F7"  # light neutral
col_high <- "#B40426"  # deep red

legend_breaks <- c(-200, -100, -50, 0, 50, 100, 200)

# Helper: nice x-axis labels using legend_labels (no underscores, rockwool/peat names)
x_lab_fun <- function(x) {
  code <- stringr::str_remove(x, "^Run_(new|old)_")
  out  <- legend_labels[code]
  out[is.na(out)] <- code[is.na(out)]
  out
}

p_hm_pub <- ggplot(hm_df_onepanel,
                   aes(x = trt_run, y = nutrient, fill = pct_capped)) +
  geom_tile() +
  geom_text(
    aes(label = ifelse(sig == "ns", "", sig)),
    size = 3.2, colour = "black", alpha = 0.75, fontface = "bold"
  ) +
  # NOTE: no geom_vline here → no dashed separator
  scale_fill_gradient2(
    name   = "∆ vs control",
    low    = col_low,
    mid    = col_mid,
    high   = col_high,
    midpoint = 0,
    limits   = c(-cap, cap),
    breaks   = legend_breaks,
    labels   = function(z) paste0(z, "%")
  ) +
  scale_x_discrete(labels = x_lab_fun) +   # <--- clean names
  labs(
    x = NULL, y = NULL,
    title    = "Nutrient change relative to each run’s control",
    subtitle = "Diverging scale centred at 0 (depletion → enrichment)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle = element_text(hjust = 0.5, colour = "grey20", margin = margin(b = 8)),
    axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y   = element_text(colour = "black"),
    panel.border  = element_rect(colour = "black", fill = NA, linewidth = 0.4),
    legend.position   = "right",
    legend.title      = element_text(face = "bold"),
    legend.key.height = unit(4.5, "mm"),
    legend.key.width  = unit(5.5, "mm"),
    plot.margin       = margin(t = 6, r = 10, b = 6, l = 6)
  )

ggsave(png_out_onepanel, p_hm_pub, width = 12, height = 7.5, dpi = 600)
ggsave(svg_out_onepanel, p_hm_pub, width = 12, height = 7.5)

print(p_hm_pub)
