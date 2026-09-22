#=============================================================
# Imputation Sensitivity Analysis
#
# Accompanies main analysis script for:
# Salivary Biomarkers for Classifying Acute Physical Fatigue: 
# A Comparative, Exploratory Study - Analysis Script
#
#
# Reproduces the manuscript's "Imputation Sensitivity Analysis" (Methods
# and Results): the primary pipeline imputes missing protein values at
# 50% of that protein's minimum detected value ("detection-limit / 2"),
# a standard left-censored/missing-not-at-random (MNAR) assumption in
# proteomics. This script tests whether the full-sample, BH-FDR-corrected
# fatigue-associated protein list is sensitive to that specific
# substitution fraction, by re-running the same screening logic under a
# sweep of three values:
#
#   1. PRIMARY  — 50% of detection limit (matches the main pipeline)
#   2. ALT-75   — 75% of detection limit
#   3. ALT-100  — 100% of detection limit (full substitution)
#=============================================================

#-------- Clear Workspace ---------
cat("\014")
rm(list = ls())

#---------- Load Libraries ----------
library(tidyverse)
library(here)
library(future)
library(furrr)
library(lme4)
library(lmerTest)

#---------- Explicit here() project root ----------
here::i_am("Scripts/PROSPER_Imputation_Sensitivity.R")

set.seed(321)

n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)
cat("Parallel backend: multisession,", n_workers, "workers\n")

#===========================================================
# SHARED SETUP — prevalence-filtered protein data
#===========================================================

cat("Loading protein data...\n")
protein_data <- read_csv(here("Data", "Fatigue_Study_10_Samples(Proteins).csv"), show_col_types = FALSE)

subject_cols <- protein_data %>%
  dplyr::select(matches("^S\\d+_T\\d+$")) %>%
  names()

df_all_proteins_long <- protein_data %>%
  dplyr::select(`Gene Symbol`, Accession, all_of(subject_cols)) %>%
  pivot_longer(
    cols      = all_of(subject_cols),
    names_to  = "Sample",
    values_to = "Abundance"
  ) %>%
  mutate(
    ID      = str_extract(Sample, "^S\\d+"),
    Measure = str_extract(Sample, "T(\\d+)$", group = 1) %>%
              as.integer() %>% as.character(),
    ID      = factor(ID),
    Measure = factor(Measure, levels = c("1","2","3","4","5","6","7"))
  ) %>%
  dplyr::select(ID, Measure, Gene_Symbol = `Gene Symbol`, Accession, Abundance)

n_samples <- 69

proteins_retained <- df_all_proteins_long %>%
  group_by(Accession) %>%
  summarise(pct_present = sum(!is.na(Abundance)) / n_samples * 100,
            .groups = "drop") %>%
  filter(pct_present > 50) %>%
  pull(Accession)

cat("Proteins retained (>50% presence):", length(proteins_retained), "\n")

harmonize_ids <- function(df) {
  df %>%
    mutate(ID = case_when(
      ID == "S01"  ~ "S001", ID == "S02"  ~ "S002",
      ID == "S03"  ~ "S003", ID == "S04"  ~ "S004",
      ID == "S05"  ~ "S005", ID == "S06"  ~ "S006",
      ID == "S07"  ~ "S007", ID == "S08"  ~ "S008",
      ID == "S09"  ~ "S009", ID == "S010" ~ "S010"
    ), ID = factor(ID)) %>%
    mutate(Measure = relevel(Measure, ref = "3"))
}

#===========================================================
# HELPER: full-sample mixed-effects screening
#===========================================================

screen_proteins_full <- function(df_imputed, protein_list) {
  future_map_dfr(protein_list, function(acc) {
    dat <- df_imputed %>% filter(Accession == acc)
    tryCatch({
      fit   <- lmer(log_Abundance ~ Measure + (1 | ID), data = dat, REML = FALSE)
      coefs <- as.data.frame(coef(summary(fit)))
      data.frame(
        Accession = acc,
        est_m4    = coefs["Measure4", "Estimate"],
        se_m4     = coefs["Measure4", "Std. Error"],
        p_fatigue = coefs["Measure4", "Pr(>|t|)"],
        p_diurnal = coefs["Measure2", "Pr(>|t|)"]
      )
    }, error = function(e) {
      data.frame(Accession = acc, est_m4 = NA_real_, se_m4 = NA_real_,
                 p_fatigue = NA_real_, p_diurnal = NA_real_)
    })
  }, .options = furrr_options(seed = TRUE, packages = c("dplyr", "lme4", "lmerTest")))
}

# Nominal fatigue-only significance, then BH-FDR correction on p_fatigue
# across the full candidate set. Returns both the nominal-only set
# (informational) and the BH-corrected set (the actual object of interest
# for this sensitivity check).
select_panel <- function(screen_df) {
  fatigue_nominal <- screen_df %>%
    filter(!is.na(p_fatigue), p_fatigue < 0.05,
           (is.na(p_diurnal) | p_diurnal > 0.05))
  screen_df_adj <- screen_df %>%
    mutate(p_fatigue_adj = p.adjust(p_fatigue, method = "BH"))
  fatigue_sig <- screen_df_adj %>%
    filter(!is.na(p_fatigue_adj), p_fatigue_adj < 0.05,
           (is.na(p_diurnal) | p_diurnal > 0.05))
  list(fatigue_nominal = fatigue_nominal, fatigue_sig = fatigue_sig,
       panel = fatigue_sig$Accession)
}

#-------- Shared helper: impute at a given fraction of detection limit --------#
impute_at_lod_fraction <- function(df_long, fraction) {
  df_long %>%
    group_by(Accession) %>%
    mutate(
      detection_limit = min(Abundance, na.rm = TRUE) * fraction,
      Abundance       = ifelse(is.na(Abundance), detection_limit, Abundance)
    ) %>%
    ungroup() %>%
    dplyr::select(-detection_limit) %>%
    mutate(log_Abundance = log(Abundance)) %>%
    harmonize_ids()
}

#===========================================================
# VARIANT 1 — PRIMARY: 50% of detection limit
#===========================================================

cat("\n=== VARIANT 1/3: PRIMARY (50% of detection limit) ===\n")

df_imputed_primary <- df_all_proteins_long %>%
  filter(Accession %in% proteins_retained) %>%
  impute_at_lod_fraction(0.50)

t0 <- Sys.time()
screen_primary <- screen_proteins_full(df_imputed_primary, unique(df_imputed_primary$Accession))
cat("  Screening time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

result_primary <- select_panel(screen_primary)
cat("  Nominal fatigue-significant proteins (p<0.05):", nrow(result_primary$fatigue_nominal), "\n")
cat("  BH-corrected fatigue-significant proteins:    ", nrow(result_primary$fatigue_sig), "\n")

#===========================================================
# VARIANT 2 — ALT-75: 75% of detection limit
#===========================================================

cat("\n=== VARIANT 2/3: ALT-75 (75% of detection limit) ===\n")

df_imputed_lod75 <- df_all_proteins_long %>%
  filter(Accession %in% proteins_retained) %>%
  impute_at_lod_fraction(0.75)

t0 <- Sys.time()
screen_lod75 <- screen_proteins_full(df_imputed_lod75, unique(df_imputed_lod75$Accession))
cat("  Screening time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

result_lod75 <- select_panel(screen_lod75)
cat("  Nominal fatigue-significant proteins (p<0.05):", nrow(result_lod75$fatigue_nominal), "\n")
cat("  BH-corrected fatigue-significant proteins:    ", nrow(result_lod75$fatigue_sig), "\n")

#===========================================================
# VARIANT 3 — ALT-100: 100% of detection limit
#===========================================================

cat("\n=== VARIANT 3/3: ALT-100 (100% of detection limit) ===\n")

df_imputed_lod100 <- df_all_proteins_long %>%
  filter(Accession %in% proteins_retained) %>%
  impute_at_lod_fraction(1.0)

t0 <- Sys.time()
screen_lod100 <- screen_proteins_full(df_imputed_lod100, unique(df_imputed_lod100$Accession))
cat("  Screening time:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")

result_lod100 <- select_panel(screen_lod100)
cat("  Nominal fatigue-significant proteins (p<0.05):", nrow(result_lod100$fatigue_nominal), "\n")
cat("  BH-corrected fatigue-significant proteins:    ", nrow(result_lod100$fatigue_sig), "\n")

#===========================================================
# COMPARISON ACROSS THE 3 IMPUTATION FRACTIONS
#===========================================================

cat("\n=== SENSITIVITY COMPARISON: 50% vs 75% vs 100% of detection limit ===\n")

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

fs_primary <- result_primary$fatigue_sig$Accession
fs_lod75   <- result_lod75$fatigue_sig$Accession
fs_lod100  <- result_lod100$fatigue_sig$Accession

cat("\nFatigue-significant protein set sizes:\n")
cat("  Primary (50% LOD):", length(fs_primary), "\n")
cat("  Alt-75  (75% LOD): ", length(fs_lod75),  "\n")
cat("  Alt-100 (100% LOD):", length(fs_lod100),  "\n")

cat("\nJaccard overlap of fatigue-significant sets:\n")
cat("  Primary vs Alt-75: ", round(jaccard(fs_primary, fs_lod75), 3), "\n")
cat("  Primary vs Alt-100:", round(jaccard(fs_primary, fs_lod100), 3), "\n")
cat("  Alt-75  vs Alt-100:", round(jaccard(fs_lod75, fs_lod100), 3), "\n")

# Gene symbol lookup for readability
gene_lookup <- protein_data %>%
  dplyr::select(Accession, `Gene Symbol`) %>%
  distinct()

# NOTE: "panel" here means the BH-corrected, full-sample fatigue-significant
# protein set (the manuscript's Imputation Sensitivity Analysis input list)
# — NOT the nested classification panels (Targeted/Protein/Combined)
# reported in Table 2 of the main analysis script.
all_bhsig_accessions <- union(union(result_primary$panel, result_lod75$panel), result_lod100$panel)

panel_comparison <- tibble(Accession = all_bhsig_accessions) %>%
  left_join(gene_lookup, by = "Accession") %>%
  mutate(
    In_Primary_BHSig = Accession %in% result_primary$panel,
    In_Alt75_BHSig   = Accession %in% result_lod75$panel,
    In_Alt100_BHSig  = Accession %in% result_lod100$panel,
    Est_m4_Primary   = result_primary$fatigue_sig$est_m4[match(Accession, result_primary$fatigue_sig$Accession)],
    Est_m4_Alt75     = result_lod75$fatigue_sig$est_m4[match(Accession, result_lod75$fatigue_sig$Accession)],
    Est_m4_Alt100    = result_lod100$fatigue_sig$est_m4[match(Accession, result_lod100$fatigue_sig$Accession)],
    N_Methods_In_BHSig = In_Primary_BHSig + In_Alt75_BHSig + In_Alt100_BHSig
  ) %>%
  arrange(desc(N_Methods_In_BHSig), Accession) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

cat("\nBH-corrected fatigue-significant protein comparison across imputation fractions:\n")
print(as.data.frame(panel_comparison), row.names = FALSE)

write_csv(panel_comparison,
          here("Results", "Imputation_Sensitivity_Panel_Comparison.csv"))

fatigue_set_summary <- tibble(
  Method                = c("Primary (50% LOD)", "Alt-75 (75% LOD)", "Alt-100 (100% LOD)"),
  N_Nominal_Significant = c(nrow(result_primary$fatigue_nominal), nrow(result_lod75$fatigue_nominal), nrow(result_lod100$fatigue_nominal)),
  N_BH_Significant      = c(length(fs_primary), length(fs_lod75), length(fs_lod100)),
  BHSig_Members         = c(paste(result_primary$panel, collapse = "; "),
                             paste(result_lod75$panel, collapse = "; "),
                             paste(result_lod100$panel, collapse = "; "))
)

write_csv(fatigue_set_summary,
          here("Results", "Imputation_Sensitivity_Summary.csv"))

cat("\nOutputs saved to Results/:\n")
cat("  Imputation_Sensitivity_Panel_Comparison.csv\n")
cat("  Imputation_Sensitivity_Summary.csv\n")
cat("\nImputation sensitivity analysis complete.\n")

plan(sequential)
