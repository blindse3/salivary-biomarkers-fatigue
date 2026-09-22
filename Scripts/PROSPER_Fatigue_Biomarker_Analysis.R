#=============================================================
# Salivary Biomarkers for Classifying Acute Physical Fatigue: 
# A Comparative, Exploratory Study - Analysis Script
# 
#
# This script reproduces the results reported in the main
# manuscript (Table 1, Table 2, Figures 2-3, permutation importance) with
# the addition of SVM-vs-ridge-logistic robustness comparison 
# required during the revision process of peer review
#
# A companion script, PROSPER_Imputation_Sensitivity.R, separately
# reproduces the imputation sensitivity analysis
#
# ---------------------------------------------------------------
# METHODOLOGICAL SUMMARY (see manuscript Methods for full detail)
# ---------------------------------------------------------------
# Panel identification (targeted, protein, and combined) uses nested
# leave-one-subject-out cross-validation (LOSOCV, k = 10): within each
# fold, candidate biomarkers are screened for a fatigue-specific response
# (Measure 3 vs. 4) independent of diurnal variation (Measure 1 vs. 2)
# using only the 9 training subjects, so no information about the held-out
# subject informs feature selection. Fatigue-response p-values are
# Benjamini-Hochberg FDR-corrected witihn each fold. A biomarker was
# included in the final reported panel only if it was retained in >= 80%
# of folds.
#
# For proteins specifically, an additional Stage-0 step narrows the 1,352
# prevalence-filtered candidates to the most reliable 50% (by rested-phase,
# Measures 1-3, measurement consistency, training subjects only) before
# any significance testing, further limiting the false-discovery burden.
#
# Classification used ridge-penalized logistic regression (glmnet,
# alpha = 0, lambda.1se via inner 5-fold CV) as the primary classifier,
# with SVM-RBF retained as a comparator from a prior version of analysis. Two
# evaluation designs are reported for every panel: Nested (each fold uses
# its own independently-selected panel which represents the unbiased, leakage-free
# estimate) and a fixed panelversion (every fold uses the same final 80%-retained
# panel, representing a secondary robustness check with a small residual leakage risk,
# since the fixed panel was chosen by looking at retention across all
# folds).
#=============================================================

#-------- Clear Workspace ---------
cat("\014")
rm(list = ls())

#---------- Load Libraries ----------
library(tidyverse)
library(here)
library(e1071)
library(pROC)
library(ggplot2)
library(ggsignif)
library(ggnewscale)
library(emmeans)
library(broom)
library(future)
library(furrr)
library(glmnet)
library(lme4)
library(lmerTest)

#---------- Explicit here() project root ----------
# Anchors here() to this script's location regardless of R's working
# directory when sourced. Adjust the relative path if this script is moved.
here::i_am("Scripts/PROSPER_Fatigue_Biomarker_Analysis.R")

#---------- Reproducibility ----------
set.seed(321)

#---------- Parallel Backend (Windows-safe) ----------
# Protein/combined-panel screening fits ~1,350+ mixed-effects models per
# fold; this parallelizes across candidates within each fold.
n_workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = n_workers)
cat("Parallel backend: multisession,", n_workers, "workers\n")

#===========================================================
# SHARED PLOT ELEMENTS
#===========================================================

phase_shading <- data.frame(
  xmin  = c(0.5, 3.5, 4.5),
  xmax  = c(3.5, 4.5, 7.5),
  phase = c("Rested", "Fatigued", "Recovery")
)

x_labels <- c("08:00", "13:00", "08:00", "12:00", "13:00", "08:00", "13:00")

day_boundaries <- c(2.5, 5.5)

#===========================================================
# SHARED HELPER FUNCTIONS
#===========================================================

# ---- Confusion-matrix -> standard classification metrics ----
# Consolidates the accuracy/sensitivity/specificity/F1/MCC/AUC computation
# used for every panel x design x classifier combination below.
compute_conf_metrics <- function(true_label, pred_label, pred_prob) {
  conf <- table(Predicted = pred_label, Actual = true_label)

  TP <- conf["Fatigued", "Fatigued"]
  TN <- conf["Rested",   "Rested"]
  FP <- conf["Fatigued", "Rested"]
  FN <- conf["Rested",   "Fatigued"]

  accuracy    <- (TP + TN) / (TP + TN + FP + FN)
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  F1          <- 2 * TP / (2 * TP + FP + FN)
  MCC         <- (TP * TN - FP * FN) /
                 sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))

  roc_obj <- roc(true_label, pred_prob, levels = c("Rested", "Fatigued"), quiet = TRUE)
  AUC     <- round(auc(roc_obj), 2)

  list(conf = conf, accuracy = accuracy, sensitivity = sensitivity,
       specificity = specificity, F1 = F1, MCC = MCC, AUC = AUC, roc = roc_obj)
}

cat_metrics <- function(label, m) {
  cat("\n", label, ":\n", sep = "")
  cat("  Accuracy:   ", round(m$accuracy, 3), "\n")
  cat("  Sensitivity:", round(m$sensitivity, 3), "\n")
  cat("  Specificity:", round(m$specificity, 3), "\n")
  cat("  F1 Score:   ", round(m$F1, 3), "\n")
  cat("  MCC:        ", round(m$MCC, 3), "\n")
  cat("  AUC:        ", m$AUC, "\n")
}

# ---- Targeted marker fold-specific selection ----
# Fits a mixed-effects model per candidate marker on training subjects
# only, extracts Tukey-adjusted p-values for the fatigue contrast
# (Measure 3 vs 4) and diurnal contrast (Measure 1 vs 2).
select_targeted_panel_fold <- function(train_ids, df_targeted_long, marker_cols) {
  purrr::map_dfr(names(marker_cols), function(m) {
    col <- marker_cols[[m]]
    dat <- df_targeted_long %>%
      filter(ID %in% train_ids) %>%
      mutate(Measure = relevel(factor(Measure), ref = "3"))

    out <- tryCatch({
      fit <- lmer(as.formula(paste0("`", col, "` ~ Measure + (1 | ID)")), data = dat)
      emm <- emmeans(fit, ~ Measure)
      pr  <- as.data.frame(pairs(emm, adjust = "tukey"))

      p_fatigue <- pr %>%
        filter(contrast %in% c("Measure3 - Measure4", "Measure4 - Measure3")) %>%
        pull(p.value) %>% .[1]
      p_diurnal <- pr %>%
        filter(contrast %in% c("Measure1 - Measure2", "Measure2 - Measure1")) %>%
        pull(p.value) %>% .[1]

      data.frame(Marker = m, p_fatigue = p_fatigue, p_diurnal = p_diurnal)
    }, error = function(e) {
      data.frame(Marker = m, p_fatigue = NA_real_, p_diurnal = NA_real_)
    })
    out
  })
}

# ---- Protein fold-specific screening (parallelized across proteins) ----
# Stage 0 (addressing alpha inflation before any significance testing):
# candidate pool is narrowed to the most reliable proteins by rested-phase
# (M1-M3) measurement consistency, computed on training subjects only.
# Default keeps the most reliable 50%.
screen_proteins_fold <- function(train_ids, df_imputed, protein_list,
                                  reliability_quantile = 0.50) {

  reliability <- df_imputed %>%
    dplyr::filter(Accession %in% protein_list, ID %in% train_ids,
                  Measure %in% c("1", "2", "3")) %>%
    dplyr::group_by(Accession, ID) %>%
    dplyr::summarise(within_subject_sd = sd(log_Abundance, na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(Accession) %>%
    dplyr::summarise(mean_rested_sd = mean(within_subject_sd, na.rm = TRUE), .groups = "drop")

  reliability_cutoff <- stats::quantile(reliability$mean_rested_sd,
                                         probs = reliability_quantile, na.rm = TRUE)

  reliable_candidates <- reliability %>%
    dplyr::filter(mean_rested_sd <= reliability_cutoff) %>%
    dplyr::pull(Accession)

  cat(sprintf("  [Stage 0: reliability pre-filter] %d of %d candidates retained (top %.0f%% by rested-phase reliability)\n",
              length(reliable_candidates), length(protein_list), reliability_quantile * 100))

  future_map_dfr(reliable_candidates, function(acc) {
    dat <- df_imputed %>% filter(Accession == acc, ID %in% train_ids)
    tryCatch({
      fit   <- lmer(log_Abundance ~ Measure + (1 | ID), data = dat, REML = FALSE)
      coefs <- as.data.frame(coef(summary(fit)))

      # p_diurnal wrapped in its own tryCatch: an emmeans/pairs() failure
      # under a near-singular fit (common with only 9 training subjects)
      # must not wipe out the already-valid est_m4/p_fatigue above it.
      p_diurnal <- tryCatch({
        emm <- emmeans::emmeans(fit, ~ Measure)
        pr  <- as.data.frame(emmeans::pairs(emm, adjust = "tukey"))
        pr %>%
          dplyr::filter(contrast %in% c("Measure1 - Measure2", "Measure2 - Measure1")) %>%
          dplyr::pull(p.value) %>% .[1]
      }, error = function(e) NA_real_)

      data.frame(
        Accession = acc,
        est_m4    = coefs["Measure4", "Estimate"],
        se_m4     = coefs["Measure4", "Std. Error"],
        p_fatigue = coefs["Measure4", "Pr(>|t|)"],
        p_diurnal = p_diurnal
      )
    }, error = function(e) {
      data.frame(Accession = acc, est_m4 = NA_real_, se_m4 = NA_real_,
                 p_fatigue = NA_real_, p_diurnal = NA_real_)
    })
  }, .options = furrr_options(seed = TRUE, packages = c("dplyr", "lme4", "lmerTest", "emmeans")))
}

# ---- Training-only z-score scaling, applied to a target set of IDs ----
build_protein_wide <- function(df_imputed, panel, scaling_ids, apply_ids) {
  scaling <- df_imputed %>%
    filter(Accession %in% panel, ID %in% scaling_ids, Measure %in% c("3", "4")) %>%
    group_by(Accession) %>%
    summarise(mean_log = mean(log_Abundance, na.rm = TRUE),
              sd_log   = sd(log_Abundance,   na.rm = TRUE),
              .groups  = "drop")

  df_imputed %>%
    filter(Accession %in% panel, ID %in% apply_ids) %>%
    left_join(scaling, by = "Accession") %>%
    mutate(z_Abundance = (log_Abundance - mean_log) / sd_log) %>%
    dplyr::select(ID, Measure, Accession, z_Abundance) %>%
    pivot_wider(names_from = Accession, values_from = z_Abundance)
}

# ---- Ridge logistic regression comparator (Reviewer C1) ----
# Fit on the SAME fold's training predictors/labels used for SVM-RBF.
fit_logistic_fold <- function(train_x, train_y, apply_x, nfolds_inner = 5) {
  x_mat <- as.matrix(train_x)
  y_num <- ifelse(as.character(train_y) == "Fatigued", 1, 0)

  # glmnet/cv.glmnet requires >= 2 predictor columns; a single-feature
  # panel (possible via the fallback path) needs ordinary logistic
  # regression instead, which handles one predictor without issue.
  if (ncol(x_mat) < 2) {
    df_train <- data.frame(y = y_num, x1 = as.numeric(x_mat[, 1]))
    fit <- tryCatch(glm(y ~ x1, data = df_train, family = "binomial"),
                     error = function(e) NULL)
    if (is.null(fit)) return(rep(NA_real_, nrow(apply_x)))
    df_apply <- data.frame(x1 = as.numeric(as.matrix(apply_x)[, 1]))
    return(as.numeric(predict(fit, newdata = df_apply, type = "response")))
  }

  cvfit <- tryCatch({
    glmnet::cv.glmnet(x_mat, y_num, family = "binomial",
                       alpha = 0, nfolds = nfolds_inner)
  }, error = function(e) NULL)

  if (is.null(cvfit)) return(rep(NA_real_, nrow(apply_x)))

  as.numeric(predict(cvfit, newx = as.matrix(apply_x),
                      s = "lambda.1se", type = "response"))
}

# ---- Retention frequency -> final reported panel (80% of folds) ----
summarize_fold_retention <- function(fold_panel_list, n_folds, threshold = 0.8) {
  all_features <- unique(unlist(fold_panel_list))
  retention <- tibble(
    Feature = all_features,
    Folds_Retained = map_int(all_features, function(f) {
      sum(map_lgl(fold_panel_list, ~ f %in% .x))
    })
  ) %>%
    mutate(Pct_Folds = Folds_Retained / n_folds) %>%
    arrange(desc(Pct_Folds))

  final_panel <- retention %>% filter(Pct_Folds >= threshold) %>% pull(Feature)
  list(retention = retention, final_panel = final_panel)
}

#=============================================================
# TABLE 1: Change in Performance Parameters
#=============================================================
# Checks that the fatiguing protocol induced acute physical fatigue: paired
# comparison of physiological/performance measures immediately before vs.
# immediately after the protocol (Measures 3 vs. 4).

cat("\n=== TABLE 1: Change in Performance Parameters ===\n")

data <- read_csv(here("Data", "Biomarker_Metadata_Saliva.csv"), show_col_types = FALSE)

performance_vars <- c(
  "Recovery_HR",
  "Queens",
  "HR",
  "HRV",
  "BR",
  "CMJ_Jump_Height(cm)",
  "IMTP_Peak_Vertical_Force(N)",
  "Max_HRPU",
  "Max_Plank(s)",
  "Anareobic_Power"
)

df_perf <- data %>%
  mutate(Group = case_when(
    Measure %in% c(1, 2, 3) ~ "Pre",
    Measure == 4            ~ "Post"
  )) %>%
  filter(!is.na(Group))

df_perf_wide <- df_perf %>%
  group_by(ID, Group) %>%
  summarise(across(all_of(performance_vars), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Group, values_from = all_of(performance_vars))

# Descriptive summary (mean/SD by Pre/Post) -> Table 1
perf_summary <- df_perf %>%
  group_by(Group) %>%
  summarise(across(all_of(performance_vars),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

cat("\nDescriptive statistics (Pre/Post):\n")
print(perf_summary %>% mutate(across(where(is.numeric), ~round(.x, 2))))

# Shapiro-Wilk normality check on each variable's paired difference score
shapiro_results <- lapply(performance_vars, function(var) {
  diff_score <- df_perf_wide[[paste0(var, "_Post")]] - df_perf_wide[[paste0(var, "_Pre")]]
  test <- tryCatch(shapiro.test(diff_score), error = function(e) NULL)
  data.frame(Variable = var,
             Shapiro_W = if (is.null(test)) NA_real_ else test$statistic,
             Shapiro_p = if (is.null(test)) NA_real_ else test$p.value)
}) %>% bind_rows()

cat("\nShapiro-Wilk normality check (paired difference scores):\n")
print(shapiro_results %>% mutate(across(where(is.numeric), ~round(.x, 3))), row.names = FALSE)

# Paired t-tests + Cohen's d / Hedges' g
perf_tests <- lapply(performance_vars, function(var) {
  pre  <- df_perf_wide[[paste0(var, "_Pre")]]
  post <- df_perf_wide[[paste0(var, "_Post")]]

  test <- t.test(pre, post, paired = TRUE)

  d_diff   <- post - pre
  n        <- sum(!is.na(d_diff))
  cohen_d  <- mean(d_diff, na.rm = TRUE) / sd(d_diff, na.rm = TRUE)
  hedges_g <- cohen_d * (1 - (3 / (4 * (n - 1) - 1)))

  broom::tidy(test) %>%
    mutate(Variable = var, n = n, Cohen_d = cohen_d, Hedges_g = hedges_g)
})

perf_results <- bind_rows(perf_tests) %>%
  dplyr::select(Variable, n, estimate, statistic, p.value, Cohen_d, Hedges_g) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

cat("\nTable 1 — Paired t-tests, Pre- vs. Post-fatigue (Measure 3 vs. 4):\n")
print(perf_results, row.names = FALSE)

write_csv(perf_results, here("Results", "Table1_Performance_Parameters.csv"))
cat("\nSaved: Results/Table1_Performance_Parameters.csv\n")

#=============================================================
# IDENTIFICATION OF CANDIDATE TARGETED BIOMARKERS
#=============================================================

cat("\n=== TARGETED BIOMARKERS: Data Preparation ===\n")

biomarkers <- c(
  "Cortisol_(ug/dL)_mean",
  "Testosterone_(pg/mL)_mean",
  "SIgA_(ug/mL)_mean",
  "Alpha-Amylase_(U/mL)_mean",
  "Uric Acid_(mg/dL)_mean",
  "IL-6_(pg/mL)"
)

df_long <- data %>%
  dplyr::select(ID, Measure, `Fatigue Status`, Time, Sex, all_of(biomarkers)) %>%
  mutate(
    ID      = factor(ID),
    Measure = factor(Measure)
  ) %>%
  mutate(
    log_cortisol     = log(`Cortisol_(ug/dL)_mean`),
    log_testosterone = log(`Testosterone_(pg/mL)_mean`),
    log_IgA          = log(`SIgA_(ug/mL)_mean`),
    log_AA           = log(`Alpha-Amylase_(U/mL)_mean`),
    log_UA           = log(`Uric Acid_(mg/dL)_mean`),
    log_IL6          = log(`IL-6_(pg/mL)`)
  )

marker_cols <- list(
  Cortisol     = "log_cortisol",
  Testosterone = "log_testosterone",
  IgA          = "log_IgA",
  AA           = "log_AA",
  UA           = "log_UA",
  IL6          = "log_IL6"
)

#-------- Temporal Analysis of Targeted Biomarker Concentrations --------#
# Full-sample mixed-effects models (all 6 candidates), Measure 3 as
# reference (most proximal resting measurement to the fatigue protocol).
# Feeds the AA/IL-6/cortisol/testosterone discussion in Results/Discussion.

df_long_ref3 <- df_long %>% mutate(Measure = relevel(Measure, ref = "3"))

targeted_models <- lapply(marker_cols, function(col) {
  lmer(as.formula(paste0("`", col, "` ~ Measure + (1 | ID)")), data = df_long_ref3)
})
names(targeted_models) <- names(marker_cols)

targeted_pairs_all <- lapply(names(targeted_models), function(m) {
  emm <- emmeans(targeted_models[[m]], ~ Measure)
  broom::tidy(pairs(emm, adjust = "tukey")) %>% mutate(Biomarker = m)
}) %>% bind_rows()

targeted_pairs_sig <- targeted_pairs_all %>%
  filter(adj.p.value < 0.05) %>%
  mutate(Pattern = case_when(
    contrast %in% c("Measure1 - Measure3", "Measure3 - Measure1") ~ "Baseline Stability",
    contrast %in% c("Measure1 - Measure2", "Measure2 - Measure1") ~ "Diurnal Effect",
    contrast %in% c("Measure3 - Measure4", "Measure4 - Measure3") ~ "Fatigue Response",
    contrast %in% c("Measure4 - Measure5", "Measure5 - Measure4") ~ "Immediate Recovery",
    contrast %in% c("Measure4 - Measure6", "Measure6 - Measure4",
                    "Measure4 - Measure7", "Measure7 - Measure4") ~ "Next Day Recovery",
    TRUE ~ "Other"
  )) %>%
  filter(Pattern != "Other") %>%
  dplyr::select(Biomarker, Pattern, contrast, estimate, std.error, statistic, adj.p.value) %>%
  arrange(Biomarker, Pattern)

cat("\nSignificant post-hoc contrasts, all 6 targeted markers:\n")
print(targeted_pairs_sig %>% mutate(across(where(is.numeric), ~round(.x, 4))), n = Inf)

write_csv(targeted_pairs_sig, here("Results", "Biomarker_Posthoc_Contrasts.csv"))
cat("Saved: Results/Biomarker_Posthoc_Contrasts.csv\n")

#-------- Part 1: Targeted Biomarker Model — Nested LOSOCV --------#
# Selection: within each fold, all 6 candidates are screened on training
# subjects only (Tukey-adjusted fatigue + diurnal contrasts, BH-FDR
# corrected within-fold); retained markers form the fold panel.

cat("\n=== IDENTIFICATION OF CANDIDATE TARGETED BIOMARKERS: Nested LOSOCV ===\n")

df_targeted_long <- df_long %>%
  dplyr::select(ID, Measure, log_cortisol, log_testosterone,
                log_IgA, log_AA, log_UA, log_IL6) %>%
  mutate(ID = factor(ID), Measure = factor(Measure, levels = c("1","2","3","4","5","6","7")))

subjects_targeted <- levels(factor(df_targeted_long$ID))
n_folds <- length(subjects_targeted)

targeted_losocv_results    <- vector("list", n_folds)
targeted_losocv_results_lr <- vector("list", n_folds)
targeted_fold_panels       <- vector("list", n_folds)
targeted_fallback_log      <- list()

for (i in seq_along(subjects_targeted)) {

  # Reseed per fold: svm(probability=TRUE)'s Platt-scaling calibration and
  # cv.glmnet's randomized fold assignment both consume random numbers, so
  # a single top-of-script set.seed() alone would not make re-running this
  # loop reproducible across sessions.
  set.seed(1000 + i)

  held_out  <- subjects_targeted[i]
  train_ids <- setdiff(subjects_targeted, held_out)

  cat(sprintf("[Targeted Nested LOSOCV] Fold %d/%d — holding out %s\n", i, n_folds, held_out))

  fold_screen <- select_targeted_panel_fold(train_ids, df_targeted_long, marker_cols) %>%
    mutate(p_fatigue_adj = p.adjust(p_fatigue, method = "BH"))

  fold_panel_names <- fold_screen %>%
    filter(!is.na(p_fatigue_adj), !is.na(p_diurnal),
           p_fatigue_adj < 0.05, p_diurnal > 0.05) %>%
    pull(Marker)

  if (length(fold_panel_names) == 0) {
    fold_panel_names <- fold_screen %>%
      filter(!is.na(p_fatigue)) %>%
      arrange(p_fatigue) %>%
      slice_head(n = 1) %>%
      pull(Marker)
    cat("  WARNING: no marker passed fatigue+diurnal criteria; falling back to smallest p_fatigue marker:",
        paste(fold_panel_names, collapse = ", "), "\n")
    targeted_fallback_log[[length(targeted_fallback_log) + 1]] <-
      data.frame(Fold = i, Held_Out = held_out, Fallback_Panel = paste(fold_panel_names, collapse = "; "))
  }

  targeted_fold_panels[[i]] <- fold_panel_names
  fold_cols <- unname(unlist(marker_cols[fold_panel_names]))
  z_cols    <- paste0("z_", fold_panel_names)

  scaling_fold <- df_targeted_long %>%
    filter(ID %in% train_ids, Measure %in% c("3", "4")) %>%
    summarise(across(all_of(fold_cols),
                      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))))

  apply_scaling <- function(dat) {
    for (j in seq_along(fold_panel_names)) {
      col  <- fold_cols[j]
      zcol <- z_cols[j]
      dat[[zcol]] <- (dat[[col]] - scaling_fold[[paste0(col, "_mean")]]) / scaling_fold[[paste0(col, "_sd")]]
    }
    dat
  }

  # Raw targeted-biomarker concentrations are never imputed (unlike
  # proteins); genuine assay-failure NAs are dropped explicitly rather than
  # silently by predict.svm(), which would otherwise desynchronize row
  # counts between true labels and predictions.
  train_m34 <- df_targeted_long %>%
    filter(ID %in% train_ids, Measure %in% c("3", "4")) %>%
    apply_scaling() %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"))) %>%
    filter(!if_any(all_of(z_cols), is.na))

  test_m34 <- df_targeted_long %>%
    filter(ID == held_out, Measure %in% c("3", "4")) %>%
    apply_scaling() %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"),
                             levels = levels(train_m34$Fatigue))) %>%
    filter(!if_any(all_of(z_cols), is.na))

  model <- svm(x = train_m34[, z_cols, drop = FALSE], y = train_m34$Fatigue,
               kernel = "radial", probability = TRUE)
  pred_m34 <- predict(model, test_m34[, z_cols, drop = FALSE], probability = TRUE)

  targeted_losocv_results[[i]] <- list(
    true_label = as.character(test_m34$Fatigue),
    pred_label = as.character(pred_m34),
    pred_probs = attr(pred_m34, "probabilities")
  )

  prob_lr_m34 <- fit_logistic_fold(train_m34[, z_cols, drop = FALSE], train_m34$Fatigue,
                                    test_m34[, z_cols, drop = FALSE])
  targeted_losocv_results_lr[[i]] <- list(
    true_label = as.character(test_m34$Fatigue),
    pred_label = ifelse(prob_lr_m34 >= 0.5, "Fatigued", "Rested"),
    pred_probs = prob_lr_m34
  )
}

cat("Targeted nested LOSOCV complete\n")

if (length(targeted_fallback_log) > 0) {
  write_csv(bind_rows(targeted_fallback_log), here("Results", "Targeted_Panel_Fallback_Log.csv"))
}

# Stability selection: >= 80% of folds,
# matching the threshold used for the protein panel below.
targeted_retention <- summarize_fold_retention(targeted_fold_panels, n_folds, threshold = 0.8)
cat("\nTargeted marker retention across folds:\n")
print(targeted_retention$retention)
write_csv(targeted_retention$retention, here("Results", "Targeted_Panel_Fold_Retention.csv"))

final_targeted_panel_names <- targeted_retention$final_panel
cat("Final targeted panel (80% of folds):", paste(final_targeted_panel_names, collapse = ", "), "\n")

targeted_true  <- unlist(lapply(targeted_losocv_results, function(x) x$true_label))
targeted_preds <- unlist(lapply(targeted_losocv_results, function(x) x$pred_label))
targeted_probs <- unlist(lapply(targeted_losocv_results, function(x) x$pred_probs[, "Fatigued"]))
targeted_metrics_svm <- compute_conf_metrics(targeted_true, targeted_preds, targeted_probs)
cat_metrics("Targeted Model Nested LOSOCV — SVM-RBF (primary)", targeted_metrics_svm)

targeted_true_lr  <- unlist(lapply(targeted_losocv_results_lr, function(x) x$true_label))
targeted_preds_lr <- unlist(lapply(targeted_losocv_results_lr, function(x) x$pred_label))
targeted_probs_lr <- unlist(lapply(targeted_losocv_results_lr, function(x) x$pred_probs))
targeted_metrics_lr <- compute_conf_metrics(targeted_true_lr, targeted_preds_lr, targeted_probs_lr)
cat_metrics("Targeted Model Nested LOSOCV — Ridge Logistic (C1 comparator)", targeted_metrics_lr)

#-------- Part 1b: Targeted Biomarker Model — Fixed-Panel LOSOCV --------#
# Secondary robustness comparator: fixes the feature set to the consensus
# panel (final_targeted_panel_names) identically across every fold, rather
# than letting each fold select its own. Carries a mild, diffuse residual
# leakage risk (the panel was chosen by looking at all-fold retention), so
# this is reported alongside, not in place of, the nested result above.

cat("\n=== IDENTIFICATION OF CANDIDATE TARGETED BIOMARKERS: Fixed-Panel LOSOCV ===\n")
cat("Fixed panel (80% of folds):", paste(final_targeted_panel_names, collapse = " + "), "\n")

fixed_targeted_cols  <- unname(unlist(marker_cols[final_targeted_panel_names]))
fixed_targeted_zcols <- paste0("z_", final_targeted_panel_names)

targeted_fixed_results    <- vector("list", n_folds)
targeted_fixed_results_lr <- vector("list", n_folds)

for (i in seq_along(subjects_targeted)) {
  set.seed(5000 + i)

  held_out  <- subjects_targeted[i]
  train_ids <- setdiff(subjects_targeted, held_out)

  scaling_fixed <- df_targeted_long %>%
    filter(ID %in% train_ids, Measure %in% c("3", "4")) %>%
    summarise(across(all_of(fixed_targeted_cols),
                      list(mean = ~mean(.x, na.rm = TRUE), sd = ~sd(.x, na.rm = TRUE))))

  apply_fixed_scaling <- function(dat) {
    for (j in seq_along(final_targeted_panel_names)) {
      col  <- fixed_targeted_cols[j]
      zcol <- fixed_targeted_zcols[j]
      dat[[zcol]] <- (dat[[col]] - scaling_fixed[[paste0(col, "_mean")]]) / scaling_fixed[[paste0(col, "_sd")]]
    }
    dat
  }

  train_fixed <- df_targeted_long %>%
    filter(ID %in% train_ids, Measure %in% c("3", "4")) %>%
    apply_fixed_scaling() %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"))) %>%
    filter(!if_any(all_of(fixed_targeted_zcols), is.na))

  test_fixed <- df_targeted_long %>%
    filter(ID == held_out, Measure %in% c("3", "4")) %>%
    apply_fixed_scaling() %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"),
                             levels = levels(train_fixed$Fatigue))) %>%
    filter(!if_any(all_of(fixed_targeted_zcols), is.na))

  model_fixed <- svm(x = train_fixed[, fixed_targeted_zcols, drop = FALSE], y = train_fixed$Fatigue,
                      kernel = "radial", probability = TRUE)
  pred_fixed <- predict(model_fixed, test_fixed[, fixed_targeted_zcols, drop = FALSE], probability = TRUE)

  targeted_fixed_results[[i]] <- list(
    true_label = as.character(test_fixed$Fatigue),
    pred_label = as.character(pred_fixed),
    pred_probs = attr(pred_fixed, "probabilities")
  )

  prob_lr_fixed <- fit_logistic_fold(train_fixed[, fixed_targeted_zcols, drop = FALSE], train_fixed$Fatigue,
                                      test_fixed[, fixed_targeted_zcols, drop = FALSE])
  targeted_fixed_results_lr[[i]] <- list(
    true_label = as.character(test_fixed$Fatigue),
    pred_label = ifelse(prob_lr_fixed >= 0.5, "Fatigued", "Rested"),
    pred_probs = prob_lr_fixed
  )
}

targeted_fx_true  <- unlist(lapply(targeted_fixed_results, function(x) x$true_label))
targeted_fx_preds <- unlist(lapply(targeted_fixed_results, function(x) x$pred_label))
targeted_fx_probs <- unlist(lapply(targeted_fixed_results, function(x) x$pred_probs[, "Fatigued"]))
targeted_metrics_svm_fx <- compute_conf_metrics(targeted_fx_true, targeted_fx_preds, targeted_fx_probs)
cat_metrics("Targeted Model Fixed-Panel LOSOCV — SVM-RBF", targeted_metrics_svm_fx)

targeted_fx_true_lr  <- unlist(lapply(targeted_fixed_results_lr, function(x) x$true_label))
targeted_fx_preds_lr <- unlist(lapply(targeted_fixed_results_lr, function(x) x$pred_label))
targeted_fx_probs_lr <- unlist(lapply(targeted_fixed_results_lr, function(x) x$pred_probs))
targeted_metrics_lr_fx <- compute_conf_metrics(targeted_fx_true_lr, targeted_fx_preds_lr, targeted_fx_probs_lr)
cat_metrics("Targeted Model Fixed-Panel LOSOCV — Ridge Logistic", targeted_metrics_lr_fx)

#=============================================================
# IDENTIFICATION OF CANDIDATE PROTEINS
#=============================================================

cat("\n=== IDENTIFICATION OF CANDIDATE PROTEINS: Data Preparation ===\n")

protein_data <- read_csv(here("Data", "Fatigue_Study_10_Samples(Proteins).csv"), show_col_types = FALSE)

subject_cols <- protein_data %>%
  dplyr::select(matches("^S\\d+_T\\d+$")) %>%
  names()

df_all_proteins_long <- protein_data %>%
  dplyr::select(`Gene Symbol`, Accession, all_of(subject_cols)) %>%
  pivot_longer(cols = all_of(subject_cols), names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    ID      = str_extract(Sample, "^S\\d+"),
    Measure = str_extract(Sample, "T(\\d+)$", group = 1) %>% as.integer() %>% as.character(),
    ID      = factor(ID),
    Measure = factor(Measure, levels = c("1","2","3","4","5","6","7"))
  ) %>%
  dplyr::select(ID, Measure, Gene_Symbol = `Gene Symbol`, Accession, Abundance)

n_samples <- 69

proteins_retained <- df_all_proteins_long %>%
  group_by(Accession) %>%
  summarise(pct_present = sum(!is.na(Abundance)) / n_samples * 100, .groups = "drop") %>%
  filter(pct_present > 50) %>%
  pull(Accession)

cat("Untargeted mass spectrometry proteins detected:", length(unique(df_all_proteins_long$Accession)), "\n")
cat("Proteins retained (>50% presence):", length(proteins_retained), "\n")

# Detection-limit imputation (50% of each protein's minimum detected
# value) + log-normalization. Label-independent, standard proteomics
# data-cleaning practice; applied once on the full sample (not nested).
df_imputed <- df_all_proteins_long %>%
  filter(Accession %in% proteins_retained) %>%
  group_by(Accession) %>%
  mutate(
    detection_limit = min(Abundance, na.rm = TRUE) * 0.5,
    Abundance       = ifelse(is.na(Abundance), detection_limit, Abundance)
  ) %>%
  ungroup() %>%
  dplyr::select(-detection_limit) %>%
  mutate(log_Abundance = log(Abundance))

#-------- Imputation Summary --------#

imputation_summary <- df_all_proteins_long %>%
  filter(Accession %in% proteins_retained) %>%
  group_by(Accession) %>%
  summarise(n_total = n(), n_missing = sum(is.na(Abundance)),
            pct_missing = round(n_missing / n_total * 100, 1), .groups = "drop")

total_values  <- nrow(df_all_proteins_long %>% filter(Accession %in% proteins_retained))
total_missing <- sum(is.na(df_all_proteins_long %>% filter(Accession %in% proteins_retained) %>% pull(Abundance)))
pct_imputed   <- round(total_missing / total_values * 100, 1)

cat("\n--- Imputation Summary ---\n")
cat("Total protein-sample values:", total_values, "\n")
cat("Missing values imputed:     ", total_missing, sprintf("(%.1f%%)\n", pct_imputed))
cat("Proteins with any missing:  ", sum(imputation_summary$n_missing > 0), "of", nrow(imputation_summary), "\n")
cat("Proteins with >10% missing: ", sum(imputation_summary$pct_missing > 10), "\n")
cat("Proteins with >20% missing: ", sum(imputation_summary$pct_missing > 20), "\n")
cat("Max missing for any protein:", max(imputation_summary$pct_missing), "%\n")

# Harmonize subject IDs to match the targeted-biomarker dataset's ID format
df_imputed <- df_imputed %>%
  mutate(ID = case_when(
    ID == "S01"  ~ "S001", ID == "S02"  ~ "S002", ID == "S03"  ~ "S003",
    ID == "S04"  ~ "S004", ID == "S05"  ~ "S005", ID == "S06"  ~ "S006",
    ID == "S07"  ~ "S007", ID == "S08"  ~ "S008", ID == "S09"  ~ "S009",
    ID == "S010" ~ "S010"
  ), ID = factor(ID)) %>%
  mutate(Measure = relevel(Measure, ref = "3"))

protein_gene_lookup <- protein_data %>%
  dplyr::select(Accession, `Gene Symbol`, Description) %>%
  distinct()

protein_list <- unique(df_imputed$Accession)

#-------- Part 2b: Protein Model — Nested LOSOCV --------#
# Selection: within each fold, Stage-0-reliable candidates are screened on
# training subjects only, BH-FDR corrected within-fold, then thresholded
# at |effect size| >= 2.0 log units to form the fold-specific panel.

cat("\n=== IDENTIFICATION OF CANDIDATE PROTEINS: Nested LOSOCV ===\n")

subjects_protein <- levels(factor(df_imputed$ID))
stopifnot(length(subjects_protein) == n_folds)

protein_losocv_results    <- vector("list", n_folds)
protein_losocv_results_lr <- vector("list", n_folds)
protein_fold_panels       <- vector("list", n_folds)
protein_fallback_log      <- list()

fold_start_time <- Sys.time()

for (i in seq_along(subjects_protein)) {
  set.seed(2000 + i)

  held_out  <- subjects_protein[i]
  train_ids <- setdiff(subjects_protein, held_out)

  cat(sprintf("\n[Protein Nested LOSOCV] Fold %d/%d — holding out %s (elapsed %.1f min)\n",
              i, n_folds, held_out, as.numeric(difftime(Sys.time(), fold_start_time, units = "mins"))))

  fold_screen <- screen_proteins_fold(train_ids, df_imputed, protein_list) %>%
    mutate(p_fatigue_adj = p.adjust(p_fatigue, method = "BH"))

  fold_fatigue <- fold_screen %>%
    filter(!is.na(p_fatigue_adj), p_fatigue_adj < 0.05,
           (is.na(p_diurnal) | p_diurnal > 0.05))

  fold_panel <- fold_fatigue %>% filter(abs(est_m4) >= 2.0) %>% pull(Accession)

  if (length(fold_panel) == 0) {
    fallback_pool <- if (nrow(fold_fatigue) > 0) fold_fatigue else fold_screen %>% filter(!is.na(est_m4))
    fold_panel <- fallback_pool %>% arrange(desc(abs(est_m4))) %>% slice_head(n = 5) %>% pull(Accession)
    cat("  WARNING: no protein met |effect size| >= 2.0 in this fold; falling back to top 5 by |effect size|\n")
    protein_fallback_log[[length(protein_fallback_log) + 1]] <-
      data.frame(Fold = i, Held_Out = held_out, Fallback_Panel = paste(fold_panel, collapse = "; "))
  }

  cat("  Fold panel (", length(fold_panel), " proteins):", paste(fold_panel, collapse = ", "), "\n")
  protein_fold_panels[[i]] <- fold_panel

  train_wide <- build_protein_wide(df_imputed, fold_panel, train_ids, train_ids) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested")))

  test_wide <- build_protein_wide(df_imputed, fold_panel, train_ids, held_out) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"),
                             levels = levels(train_wide$Fatigue)))

  model <- svm(x = train_wide[, fold_panel, drop = FALSE], y = train_wide$Fatigue,
               kernel = "radial", probability = TRUE)
  pred <- predict(model, test_wide[, fold_panel, drop = FALSE], probability = TRUE)

  protein_losocv_results[[i]] <- list(
    true_label = as.character(test_wide$Fatigue),
    pred_label = as.character(pred),
    pred_probs = attr(pred, "probabilities")
  )

  prob_lr <- fit_logistic_fold(train_wide[, fold_panel, drop = FALSE], train_wide$Fatigue,
                                test_wide[, fold_panel, drop = FALSE])
  protein_losocv_results_lr[[i]] <- list(
    true_label = as.character(test_wide$Fatigue),
    pred_label = ifelse(prob_lr >= 0.5, "Fatigued", "Rested"),
    pred_probs = prob_lr
  )
}

cat("\nProtein nested LOSOCV complete (total time:",
    round(as.numeric(difftime(Sys.time(), fold_start_time, units = "mins")), 1), "min)\n")

if (length(protein_fallback_log) > 0) {
  write_csv(bind_rows(protein_fallback_log), here("Results", "Protein_Panel_Fallback_Log.csv"))
}

protein_retention <- summarize_fold_retention(protein_fold_panels, n_folds, threshold = 0.8)
protein_retention_labeled <- protein_retention$retention %>%
  left_join(protein_gene_lookup, by = c("Feature" = "Accession")) %>%
  arrange(desc(Pct_Folds))

cat("\nProtein retention across folds (top 20):\n")
print(head(protein_retention_labeled, 20))
write_csv(protein_retention_labeled, here("Results", "Protein_Panel_Fold_Retention.csv"))

final_protein_panel <- protein_retention$final_panel
cat("Final protein panel (80% of folds):", paste(final_protein_panel, collapse = ", "), "\n")

protein_true  <- unlist(lapply(protein_losocv_results, function(x) x$true_label))
protein_preds <- unlist(lapply(protein_losocv_results, function(x) x$pred_label))
protein_probs <- unlist(lapply(protein_losocv_results, function(x) x$pred_probs[, "Fatigued"]))
protein_metrics_svm <- compute_conf_metrics(protein_true, protein_preds, protein_probs)
cat_metrics("Protein Model Nested LOSOCV — SVM-RBF (primary)", protein_metrics_svm)

protein_true_lr  <- unlist(lapply(protein_losocv_results_lr, function(x) x$true_label))
protein_preds_lr <- unlist(lapply(protein_losocv_results_lr, function(x) x$pred_label))
protein_probs_lr <- unlist(lapply(protein_losocv_results_lr, function(x) x$pred_probs))
protein_metrics_lr <- compute_conf_metrics(protein_true_lr, protein_preds_lr, protein_probs_lr)
cat_metrics("Protein Model Nested LOSOCV — Ridge Logistic (C1 comparator)", protein_metrics_lr)

#-------- Part 2d: Protein Model — Fixed-Panel LOSOCV --------#

cat("\n=== IDENTIFICATION OF CANDIDATE PROTEINS: Fixed-Panel LOSOCV ===\n")
cat("Fixed panel (80% of folds):", paste(final_protein_panel, collapse = ", "), "\n")

protein_fixed_results    <- vector("list", n_folds)
protein_fixed_results_lr <- vector("list", n_folds)

for (i in seq_along(subjects_protein)) {
  set.seed(6000 + i)

  held_out  <- subjects_protein[i]
  train_ids <- setdiff(subjects_protein, held_out)

  train_wide_fixed <- build_protein_wide(df_imputed, final_protein_panel, train_ids, train_ids) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested")))

  test_wide_fixed <- build_protein_wide(df_imputed, final_protein_panel, train_ids, held_out) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"),
                             levels = levels(train_wide_fixed$Fatigue)))

  model_fixed <- svm(x = train_wide_fixed[, final_protein_panel, drop = FALSE], y = train_wide_fixed$Fatigue,
                      kernel = "radial", probability = TRUE)
  pred_fixed <- predict(model_fixed, test_wide_fixed[, final_protein_panel, drop = FALSE], probability = TRUE)

  protein_fixed_results[[i]] <- list(
    true_label = as.character(test_wide_fixed$Fatigue),
    pred_label = as.character(pred_fixed),
    pred_probs = attr(pred_fixed, "probabilities")
  )

  prob_lr_fixed <- fit_logistic_fold(train_wide_fixed[, final_protein_panel, drop = FALSE], train_wide_fixed$Fatigue,
                                      test_wide_fixed[, final_protein_panel, drop = FALSE])
  protein_fixed_results_lr[[i]] <- list(
    true_label = as.character(test_wide_fixed$Fatigue),
    pred_label = ifelse(prob_lr_fixed >= 0.5, "Fatigued", "Rested"),
    pred_probs = prob_lr_fixed
  )
}

protein_fx_true  <- unlist(lapply(protein_fixed_results, function(x) x$true_label))
protein_fx_preds <- unlist(lapply(protein_fixed_results, function(x) x$pred_label))
protein_fx_probs <- unlist(lapply(protein_fixed_results, function(x) x$pred_probs[, "Fatigued"]))
protein_metrics_svm_fx <- compute_conf_metrics(protein_fx_true, protein_fx_preds, protein_fx_probs)
cat_metrics("Protein Model Fixed-Panel LOSOCV — SVM-RBF", protein_metrics_svm_fx)

protein_fx_true_lr  <- unlist(lapply(protein_fixed_results_lr, function(x) x$true_label))
protein_fx_preds_lr <- unlist(lapply(protein_fixed_results_lr, function(x) x$pred_label))
protein_fx_probs_lr <- unlist(lapply(protein_fixed_results_lr, function(x) x$pred_probs))
protein_metrics_lr_fx <- compute_conf_metrics(protein_fx_true_lr, protein_fx_preds_lr, protein_fx_probs_lr)
cat_metrics("Protein Model Fixed-Panel LOSOCV — Ridge Logistic", protein_metrics_lr_fx)

#=============================================================
# IDENTIFICATION OF COMBINED CANDIDATE BIOMARKER PANEL
#=============================================================
# Targeted markers and proteins are screened jointly within the same
# nested LOSOCV framework, rather than combining each platform's
# separately-identified panel post hoc. Because effect-size/significance
# estimates differ systematically in scale between assay types, candidates
# are ranked by a scale-invariant t-statistic and capped at the 5
# highest-ranked per fold, rather than applying one magnitude threshold
# across both assay types.

cat("\n=== IDENTIFICATION OF COMBINED PANEL: Nested LOSOCV (Unified Joint Screen) ===\n")

df_targeted_as_schema <- df_targeted_long %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, log_cortisol, log_testosterone, log_IgA, log_AA, log_UA, log_IL6) %>%
  pivot_longer(cols = starts_with("log_"), names_to = "Accession", values_to = "log_Abundance") %>%
  mutate(Accession = paste0("TARGETED_", sub("^log_", "", Accession)))

df_protein_as_schema <- df_imputed %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, Accession, log_Abundance)

df_combined_candidates <- bind_rows(df_targeted_as_schema, df_protein_as_schema) %>%
  mutate(ID = factor(ID), Measure = factor(Measure, levels = c("1","2","3","4","5","6","7"))) %>%
  mutate(Measure = relevel(Measure, ref = "3"))

# Maps a targeted marker's short name (e.g. "IgA") to its schema Accession
# (e.g. "TARGETED_IgA"), for cross-referencing fold panels against the
# unified candidate pool's naming convention.
targeted_name_to_schema <- setNames(
  paste0("TARGETED_", sub("^log_", "", unlist(marker_cols))),
  names(marker_cols)
)

combined_candidate_list <- unique(df_combined_candidates$Accession)
cat("Combined candidate pool:", length(combined_candidate_list), "features (",
    sum(grepl("^TARGETED_", combined_candidate_list)), "targeted +",
    sum(!grepl("^TARGETED_", combined_candidate_list)), "proteins)\n")

subjects_unified <- levels(factor(df_combined_candidates$ID))
stopifnot(length(subjects_unified) == n_folds)

unified_losocv_results    <- vector("list", n_folds)
unified_losocv_results_lr <- vector("list", n_folds)
unified_fold_panels       <- vector("list", n_folds)
unified_fallback_log      <- list()

fold_start_time_unified <- Sys.time()

for (i in seq_along(subjects_unified)) {
  set.seed(10000 + i)

  held_out  <- subjects_unified[i]
  train_ids <- setdiff(subjects_unified, held_out)

  cat(sprintf("\n[Combined Nested LOSOCV] Fold %d/%d — holding out %s (elapsed %.1f min)\n",
              i, n_folds, held_out, as.numeric(difftime(Sys.time(), fold_start_time_unified, units = "mins"))))

  fold_screen <- screen_proteins_fold(train_ids, df_combined_candidates, combined_candidate_list) %>%
    mutate(t_stat = est_m4 / se_m4, p_fatigue_adj = p.adjust(p_fatigue, method = "BH"))

  fold_fatigue <- fold_screen %>%
    filter(!is.na(p_fatigue_adj), p_fatigue_adj < 0.05,
           (is.na(p_diurnal) | p_diurnal > 0.05))

  # Rank-and-cap (top 5 by |t-stat|), not a fixed threshold: |t| rewards
  # measurement precision, not effect magnitude, and a fixed cutoff can
  # admit far too many candidates when many proteins clear it on precision
  # alone. Capping panel SIZE directly avoids that failure mode.
  fold_panel <- fold_fatigue %>% arrange(desc(abs(t_stat))) %>% slice_head(n = 5) %>% pull(Accession)

  if (length(fold_panel) == 0) {
    fold_panel <- fold_screen %>% filter(!is.na(t_stat)) %>%
      arrange(desc(abs(t_stat))) %>% slice_head(n = 5) %>% pull(Accession)
    cat("  WARNING: no feature was nominally significant in this fold; falling back to top 5 by |t| regardless of significance\n")
    unified_fallback_log[[length(unified_fallback_log) + 1]] <-
      data.frame(Fold = i, Held_Out = held_out, Fallback_Panel = paste(fold_panel, collapse = "; "))
  }

  cat("  Fold panel (", length(fold_panel), " features):", paste(fold_panel, collapse = ", "), "\n")
  unified_fold_panels[[i]] <- fold_panel

  # Targeted markers (unlike imputed proteins) can have genuine missing
  # data; drop NA rows explicitly before training/predicting.
  train_wide <- build_protein_wide(df_combined_candidates, fold_panel, train_ids, train_ids) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"))) %>%
    filter(!if_any(all_of(fold_panel), is.na))

  test_wide <- build_protein_wide(df_combined_candidates, fold_panel, train_ids, held_out) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"),
                             levels = levels(train_wide$Fatigue))) %>%
    filter(!if_any(all_of(fold_panel), is.na))

  model <- svm(x = train_wide[, fold_panel, drop = FALSE], y = train_wide$Fatigue,
               kernel = "radial", probability = TRUE)
  pred <- predict(model, test_wide[, fold_panel, drop = FALSE], probability = TRUE)

  unified_losocv_results[[i]] <- list(
    true_label = as.character(test_wide$Fatigue),
    pred_label = as.character(pred),
    pred_probs = attr(pred, "probabilities")
  )

  prob_lr <- fit_logistic_fold(train_wide[, fold_panel, drop = FALSE], train_wide$Fatigue,
                                test_wide[, fold_panel, drop = FALSE])
  unified_losocv_results_lr[[i]] <- list(
    true_label = as.character(test_wide$Fatigue),
    pred_label = ifelse(prob_lr >= 0.5, "Fatigued", "Rested"),
    pred_probs = prob_lr
  )
}

cat("\nCombined nested LOSOCV complete (total time:",
    round(as.numeric(difftime(Sys.time(), fold_start_time_unified, units = "mins")), 1), "min)\n")

if (length(unified_fallback_log) > 0) {
  write_csv(bind_rows(unified_fallback_log), here("Results", "Combined_Panel_Fallback_Log.csv"))
}

unified_retention <- summarize_fold_retention(unified_fold_panels, n_folds, threshold = 0.8)
unified_retention_labeled <- unified_retention$retention %>%
  left_join(protein_gene_lookup, by = c("Feature" = "Accession")) %>%
  arrange(desc(Pct_Folds))

cat("\nCombined-panel candidate retention across folds (top 20):\n")
print(head(unified_retention_labeled, 20))
write_csv(unified_retention_labeled, here("Results", "Combined_Panel_Fold_Retention.csv"))

final_combined_panel <- unified_retention$final_panel
cat("\nFinal combined panel (80% of folds):", paste(final_combined_panel, collapse = ", "), "\n")

unified_true  <- unlist(lapply(unified_losocv_results, function(x) x$true_label))
unified_preds <- unlist(lapply(unified_losocv_results, function(x) x$pred_label))
unified_probs <- unlist(lapply(unified_losocv_results, function(x) x$pred_probs[, "Fatigued"]))
combined_metrics_svm <- compute_conf_metrics(unified_true, unified_preds, unified_probs)
cat_metrics("Combined Model Nested LOSOCV — SVM-RBF (primary)", combined_metrics_svm)

unified_true_lr  <- unlist(lapply(unified_losocv_results_lr, function(x) x$true_label))
unified_preds_lr <- unlist(lapply(unified_losocv_results_lr, function(x) x$pred_label))
unified_probs_lr <- unlist(lapply(unified_losocv_results_lr, function(x) x$pred_probs))
combined_metrics_lr <- compute_conf_metrics(unified_true_lr, unified_preds_lr, unified_probs_lr)
cat_metrics("Combined Model Nested LOSOCV — Ridge Logistic (C1 comparator)", combined_metrics_lr)

#-------- Combined Panel — Fixed-Panel LOSOCV --------#

cat("\n=== IDENTIFICATION OF COMBINED PANEL: Fixed-Panel LOSOCV ===\n")
cat("Fixed panel (80% of folds):", paste(final_combined_panel, collapse = " + "), "\n")

combined_fixed_results    <- vector("list", n_folds)
combined_fixed_results_lr <- vector("list", n_folds)

for (i in seq_along(subjects_unified)) {
  set.seed(11000 + i)

  held_out  <- subjects_unified[i]
  train_ids <- setdiff(subjects_unified, held_out)

  train_wide_fixed <- build_protein_wide(df_combined_candidates, final_combined_panel, train_ids, train_ids) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested")))

  test_wide_fixed <- build_protein_wide(df_combined_candidates, final_combined_panel, train_ids, held_out) %>%
    filter(Measure %in% c("3", "4")) %>%
    mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"),
                             levels = levels(train_wide_fixed$Fatigue)))

  model_fixed <- svm(x = train_wide_fixed[, final_combined_panel, drop = FALSE], y = train_wide_fixed$Fatigue,
                      kernel = "radial", probability = TRUE)
  pred_fixed <- predict(model_fixed, test_wide_fixed[, final_combined_panel, drop = FALSE], probability = TRUE)

  combined_fixed_results[[i]] <- list(
    true_label = as.character(test_wide_fixed$Fatigue),
    pred_label = as.character(pred_fixed),
    pred_probs = attr(pred_fixed, "probabilities")
  )

  prob_lr_fixed <- fit_logistic_fold(train_wide_fixed[, final_combined_panel, drop = FALSE], train_wide_fixed$Fatigue,
                                      test_wide_fixed[, final_combined_panel, drop = FALSE])
  combined_fixed_results_lr[[i]] <- list(
    true_label = as.character(test_wide_fixed$Fatigue),
    pred_label = ifelse(prob_lr_fixed >= 0.5, "Fatigued", "Rested"),
    pred_probs = prob_lr_fixed
  )
}

combined_fx_true_svm  <- unlist(lapply(combined_fixed_results, function(x) x$true_label))
combined_fx_pred_svm  <- unlist(lapply(combined_fixed_results, function(x) x$pred_label))
combined_fx_probs_svm <- unlist(lapply(combined_fixed_results, function(x) x$pred_probs[, "Fatigued"]))
combined_metrics_svm_fx <- compute_conf_metrics(combined_fx_true_svm, combined_fx_pred_svm, combined_fx_probs_svm)
cat_metrics("Combined Model Fixed-Panel LOSOCV — SVM-RBF", combined_metrics_svm_fx)

combined_fx_true_lr  <- unlist(lapply(combined_fixed_results_lr, function(x) x$true_label))
combined_fx_pred_lr  <- unlist(lapply(combined_fixed_results_lr, function(x) x$pred_label))
combined_fx_probs_lr <- unlist(lapply(combined_fixed_results_lr, function(x) x$pred_probs))
combined_metrics_lr_fx <- compute_conf_metrics(combined_fx_true_lr, combined_fx_pred_lr, combined_fx_probs_lr)
cat_metrics("Combined Model Fixed-Panel LOSOCV — Ridge Logistic", combined_metrics_lr_fx)

#=============================================================
# TABLE 2: Model Performance Comparison
#=============================================================
# Consolidates the ridge-logistic-regression nested and fixed-panel
# results for all three panels (targeted, protein, combined) into the
# manuscript's Table 2. Sensitivity/specificity/accuracy are reported as
# both percentages and raw counts (every design evaluates a balanced
# 10 rested vs. 10 fatigued set).

cat("\n=== TABLE 2: Model Performance Comparison ===\n")

fmt_pct_count <- function(pct, denom) sprintf("%.0f%% (%d/%d)", pct * 100, round(pct * denom), denom)

table2_performance <- tibble::tribble(
  ~Panel,     ~Design,       ~Sensitivity,                    ~Specificity,                    ~Accuracy,                    ~F1,                              ~MCC,                              ~AUC,
  "Targeted", "Nested",      targeted_metrics_lr$sensitivity,    targeted_metrics_lr$specificity,    targeted_metrics_lr$accuracy,    targeted_metrics_lr$F1,             targeted_metrics_lr$MCC,             targeted_metrics_lr$AUC,
  "Targeted", "Fixed Panel", targeted_metrics_lr_fx$sensitivity, targeted_metrics_lr_fx$specificity, targeted_metrics_lr_fx$accuracy, targeted_metrics_lr_fx$F1,          targeted_metrics_lr_fx$MCC,          targeted_metrics_lr_fx$AUC,
  "Protein",  "Nested",      protein_metrics_lr$sensitivity,     protein_metrics_lr$specificity,     protein_metrics_lr$accuracy,     protein_metrics_lr$F1,              protein_metrics_lr$MCC,              protein_metrics_lr$AUC,
  "Protein",  "Fixed Panel", protein_metrics_lr_fx$sensitivity,  protein_metrics_lr_fx$specificity,  protein_metrics_lr_fx$accuracy,  protein_metrics_lr_fx$F1,           protein_metrics_lr_fx$MCC,           protein_metrics_lr_fx$AUC,
  "Combined", "Nested",      combined_metrics_lr$sensitivity,    combined_metrics_lr$specificity,    combined_metrics_lr$accuracy,    combined_metrics_lr$F1,             combined_metrics_lr$MCC,             combined_metrics_lr$AUC,
  "Combined", "Fixed Panel", combined_metrics_lr_fx$sensitivity, combined_metrics_lr_fx$specificity, combined_metrics_lr_fx$accuracy, combined_metrics_lr_fx$F1,          combined_metrics_lr_fx$MCC,          combined_metrics_lr_fx$AUC
) %>%
  mutate(
    Sensitivity_fmt = fmt_pct_count(Sensitivity, 10),
    Specificity_fmt = fmt_pct_count(Specificity, 10),
    Accuracy_fmt    = fmt_pct_count(Accuracy, 20),
    F1  = round(F1, 3), MCC = round(MCC, 3), AUC = round(AUC, 3)
  ) %>%
  dplyr::select(Panel, Design, Sensitivity = Sensitivity_fmt, Specificity = Specificity_fmt,
                Accuracy = Accuracy_fmt, F1, MCC, AUC)

cat("\nTable 2 — Model Performance Comparison (Ridge Logistic Regression):\n")
print(as.data.frame(table2_performance), row.names = FALSE)

write_csv(table2_performance, here("Results", "Table2_Model_Performance_Comparison.csv"))
cat("Saved: Results/Table2_Model_Performance_Comparison.csv\n")

#=============================================================
# SVM-RBF vs. RIDGE LOGISTIC REGRESSION COMPARISON
#=============================================================
# SVM-RBF was the originally-submitted classifier. In response to a
# reviewer concern about model flexibility relative to sample size, a
# simpler ridge-penalized logistic regression was added as a comparator,
# evaluated on the identical nested LOSOCV folds and fold-specific panels.
# Ridge logistic regression was adopted as the primary/reported classifier
# throughout the manuscript (see Table 2); this comparison is retained as
# a supplementary robustness check.

cat("\n=== SVM-RBF vs. Ridge Logistic Regression (Response to Reviewer, C1) ===\n")

model_comparison <- tibble::tribble(
  ~Panel,      ~Classifier,        ~Accuracy,                        ~Sensitivity,                        ~Specificity,                        ~F1,                        ~MCC,                        ~AUC,
  "Targeted",  "SVM-RBF",          targeted_metrics_svm$accuracy,    targeted_metrics_svm$sensitivity,    targeted_metrics_svm$specificity,    targeted_metrics_svm$F1,    targeted_metrics_svm$MCC,    targeted_metrics_svm$AUC,
  "Targeted",  "Ridge Logistic",   targeted_metrics_lr$accuracy,     targeted_metrics_lr$sensitivity,     targeted_metrics_lr$specificity,     targeted_metrics_lr$F1,     targeted_metrics_lr$MCC,     targeted_metrics_lr$AUC,
  "Protein",   "SVM-RBF",          protein_metrics_svm$accuracy,     protein_metrics_svm$sensitivity,     protein_metrics_svm$specificity,     protein_metrics_svm$F1,     protein_metrics_svm$MCC,     protein_metrics_svm$AUC,
  "Protein",   "Ridge Logistic",   protein_metrics_lr$accuracy,      protein_metrics_lr$sensitivity,      protein_metrics_lr$specificity,      protein_metrics_lr$F1,      protein_metrics_lr$MCC,      protein_metrics_lr$AUC,
  "Combined",  "SVM-RBF",          combined_metrics_svm$accuracy,    combined_metrics_svm$sensitivity,    combined_metrics_svm$specificity,    combined_metrics_svm$F1,    combined_metrics_svm$MCC,    combined_metrics_svm$AUC,
  "Combined",  "Ridge Logistic",   combined_metrics_lr$accuracy,     combined_metrics_lr$sensitivity,     combined_metrics_lr$specificity,     combined_metrics_lr$F1,     combined_metrics_lr$MCC,     combined_metrics_lr$AUC
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(as.data.frame(model_comparison), row.names = FALSE)
write_csv(model_comparison, here("Results", "Model_Comparison_SVM_vs_Logistic.csv"))

# Supplementary: fixed-panel single-protein (C1S) SVM vs. Logistic, since
# the protein panel's fixed design IS the single-protein (C1S) panel.
diagnostic_c1s_comparison <- tibble::tribble(
  ~Panel,      ~Classifier,       ~Accuracy,                          ~Sensitivity,                          ~Specificity,                          ~F1,                          ~MCC,                          ~AUC,
  "C1S only",  "SVM-RBF",         protein_metrics_svm_fx$accuracy,    protein_metrics_svm_fx$sensitivity,    protein_metrics_svm_fx$specificity,    protein_metrics_svm_fx$F1,    protein_metrics_svm_fx$MCC,    protein_metrics_svm_fx$AUC,
  "C1S only",  "Ridge Logistic",  protein_metrics_lr_fx$accuracy,     protein_metrics_lr_fx$sensitivity,     protein_metrics_lr_fx$specificity,     protein_metrics_lr_fx$F1,     protein_metrics_lr_fx$MCC,     protein_metrics_lr_fx$AUC
) %>%
  mutate(across(where(is.numeric), ~round(.x, 3)))

print(as.data.frame(diagnostic_c1s_comparison), row.names = FALSE)
write_csv(diagnostic_c1s_comparison, here("Results", "Diagnostic_C1SOnly_FixedPanel_Comparison.csv"))

#-------- Supplementary ROC figures --------#

roc_df_comparison <- bind_rows(
  data.frame(FPR = 1 - targeted_metrics_svm$roc$specificities, TPR = targeted_metrics_svm$roc$sensitivities,
             Panel = "Targeted", Classifier = "SVM-RBF") %>% arrange(FPR, TPR),
  data.frame(FPR = 1 - targeted_metrics_lr$roc$specificities,  TPR = targeted_metrics_lr$roc$sensitivities,
             Panel = "Targeted", Classifier = "Ridge Logistic") %>% arrange(FPR, TPR),
  data.frame(FPR = 1 - protein_metrics_svm$roc$specificities,  TPR = protein_metrics_svm$roc$sensitivities,
             Panel = "Protein",  Classifier = "SVM-RBF") %>% arrange(FPR, TPR),
  data.frame(FPR = 1 - protein_metrics_lr$roc$specificities,   TPR = protein_metrics_lr$roc$sensitivities,
             Panel = "Protein",  Classifier = "Ridge Logistic") %>% arrange(FPR, TPR),
  data.frame(FPR = 1 - combined_metrics_svm$roc$specificities, TPR = combined_metrics_svm$roc$sensitivities,
             Panel = "Combined", Classifier = "SVM-RBF") %>% arrange(FPR, TPR),
  data.frame(FPR = 1 - combined_metrics_lr$roc$specificities,  TPR = combined_metrics_lr$roc$sensitivities,
             Panel = "Combined", Classifier = "Ridge Logistic") %>% arrange(FPR, TPR)
)

p_roc_comparison <- ggplot(roc_df_comparison, aes(x = FPR, y = TPR, color = Panel, linetype = Classifier)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.8) +
  geom_step(linewidth = 1.1, direction = "hv") +
  scale_color_manual(values = c("Targeted" = "#E69F00", "Protein" = "#0072B2", "Combined" = "#009E73")) +
  scale_linetype_manual(values = c("SVM-RBF" = "solid", "Ridge Logistic" = "dashed")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.05, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(x = "False Positive Rate", y = "True Positive Rate", color = "Panel", linetype = "Classifier",
       title = "SVM-RBF vs. Ridge Logistic Regression (Nested LOSOCV)") +
  theme_classic() +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5), legend.position = "bottom")

ggsave(here("Results", "Model_Comparison_ROC.tiff"), plot = p_roc_comparison,
       width = 7, height = 7, units = "in", dpi = 300)
cat("Saved: Results/Model_Comparison_ROC.tiff (supplementary)\n")

#=============================================================
# PERMUTATION IMPORTANCE: Combined Model
#=============================================================
# Assesses each biomarker's relative contribution to the combined model's
# classification accuracy. A ridge-penalized logistic regression is fit to
# the full sample (Measures 3-4); for each biomarker, values are randomly
# permuted across observations (1,000 iterations) and the fitted model is
# reapplied without refitting. Importance = mean reduction in accuracy
# relative to the unpermuted baseline.

cat("\n=== PERMUTATION IMPORTANCE: Combined Model (full sample) ===\n")

# C1S accession (NP_958850.1) is hardcoded here, matching its use in
# Figure 2 below — both draw on the same specific protein, independent of
# whatever else does or doesn't end up in final_protein_panel.
df_c1s_long <- df_imputed %>%
  filter(Accession == "NP_958850.1") %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, log_C1S = log_Abundance)

df_ang_long <- df_imputed %>%
  filter(Accession == "NP_001091046.1") %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, log_ANG = log_Abundance)

df_combined_full <- df_targeted_long %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, log_IgA, log_UA) %>%
  inner_join(df_c1s_long, by = c("ID", "Measure")) %>%
  inner_join(df_ang_long, by = c("ID", "Measure")) %>%
  filter(Measure %in% c("3", "4")) %>%
  mutate(
    z_IgA   = as.numeric(scale(log_IgA)),
    z_UA    = as.numeric(scale(log_UA)),
    z_C1S   = as.numeric(scale(log_C1S)),
    z_ANG   = as.numeric(scale(log_ANG)),
    Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"))
  )

combined_panel_cols <- c("z_IgA", "z_UA", "z_C1S", "z_ANG")

set.seed(4343)
combined_logistic_full <- glmnet::cv.glmnet(
  x = as.matrix(df_combined_full[, combined_panel_cols, drop = FALSE]),
  y = ifelse(df_combined_full$Fatigue == "Fatigued", 1, 0),
  family = "binomial", alpha = 0, nfolds = 5
)

predict_combined_class <- function(newdata) {
  probs <- as.numeric(predict(combined_logistic_full,
                               newx = as.matrix(newdata[, combined_panel_cols, drop = FALSE]),
                               s = "lambda.1se", type = "response"))
  factor(ifelse(probs >= 0.5, "Fatigued", "Rested"), levels = levels(df_combined_full$Fatigue))
}

baseline_acc_combined <- mean(predict_combined_class(df_combined_full) == df_combined_full$Fatigue)

n_perm <- 1000
set.seed(123)
importance_combined <- sapply(combined_panel_cols, function(col) {
  perm_accs <- replicate(n_perm, {
    df_perm <- df_combined_full
    df_perm[[col]] <- sample(df_perm[[col]])
    mean(predict_combined_class(df_perm) == df_perm$Fatigue)
  })
  baseline_acc_combined - mean(perm_accs)
})

combined_marker_labels <- c(z_IgA = "IgA", z_UA = "UA", z_C1S = "C1S", z_ANG = "ANG")

importance_combined_df <- data.frame(
  Marker     = unname(combined_marker_labels[combined_panel_cols]),
  Importance = round(importance_combined, 3)
) %>% arrange(desc(Importance))

cat("\nCombined Model Permutation Importance (full sample, panel: IgA + UA + C1S + ANG):\n")
print(importance_combined_df, row.names = FALSE)

write_csv(importance_combined_df, here("Results", "Combined_Model_Permutation_Importance.csv"))
cat("Saved: Results/Combined_Model_Permutation_Importance.csv\n")

#=============================================================
# PREDICTED-FATIGUED CLASSIFICATION RATE TRAJECTORIES (Figure 3)
#=============================================================
# Applies each panel's ridge-logistic classifier across all 7 collection
# timepoints for every held-out subject, to assess temporal generalizability
# beyond the Measure 3/4 training window. NESTED: each held-out subject's
# probabilities come from the fold where they were excluded, using that
# fold's independently-selected panel. FIXED PANEL: same LOSOCV structure,
# but the same fixed panel identity is used in every fold.

cat("\n=== PREDICTED-FATIGUED CLASSIFICATION RATE TRAJECTORIES (Figure 3) ===\n")

compute_trajectory_losocv <- function(subjects, panel_per_fold, df_source, seed_base) {
  purrr::map_dfr(seq_along(subjects), function(i) {
    set.seed(seed_base + i)

    held_out  <- subjects[i]
    train_ids <- setdiff(subjects, held_out)
    panel     <- panel_per_fold[[i]]

    train_wide <- build_protein_wide(df_source, panel, train_ids, train_ids) %>%
      filter(Measure %in% c("3", "4")) %>%
      mutate(Fatigue = factor(ifelse(Measure == "4", "Fatigued", "Rested"))) %>%
      filter(!if_any(all_of(panel), is.na))

    test_wide_all <- build_protein_wide(df_source, panel, train_ids, held_out) %>%
      filter(!if_any(all_of(panel), is.na))

    if (nrow(train_wide) < 2 || nrow(test_wide_all) == 0) return(NULL)

    prob <- fit_logistic_fold(train_wide[, panel, drop = FALSE], train_wide$Fatigue,
                               test_wide_all[, panel, drop = FALSE])

    test_wide_all %>% mutate(Fatigue_Prob = prob) %>% dplyr::select(ID, Measure, Fatigue_Prob)
  })
}

targeted_panel_per_fold_nested <- lapply(targeted_fold_panels, function(p) unname(targeted_name_to_schema[p]))
protein_panel_per_fold_nested  <- protein_fold_panels
combined_panel_per_fold_nested <- unified_fold_panels

targeted_fixed_panel <- unname(targeted_name_to_schema[final_targeted_panel_names])
protein_fixed_panel  <- final_protein_panel
combined_fixed_panel <- final_combined_panel

targeted_panel_per_fold_fixed <- rep(list(targeted_fixed_panel), n_folds)
protein_panel_per_fold_fixed  <- rep(list(protein_fixed_panel), n_folds)
combined_panel_per_fold_fixed <- rep(list(combined_fixed_panel), n_folds)

traj_targeted_nested <- compute_trajectory_losocv(subjects_targeted, targeted_panel_per_fold_nested, df_combined_candidates, seed_base = 20001)
traj_targeted_fixed  <- compute_trajectory_losocv(subjects_targeted, targeted_panel_per_fold_fixed,  df_combined_candidates, seed_base = 20101)
traj_protein_nested  <- compute_trajectory_losocv(subjects_protein, protein_panel_per_fold_nested, df_combined_candidates, seed_base = 20201)
traj_protein_fixed   <- compute_trajectory_losocv(subjects_protein, protein_panel_per_fold_fixed,  df_combined_candidates, seed_base = 20301)
traj_combined_nested <- compute_trajectory_losocv(subjects_unified, combined_panel_per_fold_nested, df_combined_candidates, seed_base = 20401)
traj_combined_fixed  <- compute_trajectory_losocv(subjects_unified, combined_panel_per_fold_fixed,  df_combined_candidates, seed_base = 20501)

summarize_traj_accuracy <- function(df, panel_name, design_name) {
  df %>%
    mutate(Measure = factor(Measure, levels = c("1","2","3","4","5","6","7")),
           predicted_fatigued = Fatigue_Prob >= 0.5) %>%
    group_by(Measure) %>%
    summarise(pct_correct = mean(predicted_fatigued, na.rm = TRUE),
              n_subjects  = sum(!is.na(predicted_fatigued)),
              n_fatigued  = sum(predicted_fatigued, na.rm = TRUE),
              .groups = "drop") %>%
    # Wilson-score 95% CI (prop.test): this is a binomial proportion out of
    # n_subjects, not a continuous measurement, so an SD would just be a
    # fixed function of the point estimate rather than new information.
    rowwise() %>%
    mutate(
      ci_low  = if (n_subjects > 0) suppressWarnings(prop.test(n_fatigued, n_subjects)$conf.int[1]) else NA_real_,
      ci_high = if (n_subjects > 0) suppressWarnings(prop.test(n_fatigued, n_subjects)$conf.int[2]) else NA_real_
    ) %>%
    ungroup() %>%
    mutate(Panel = panel_name, Design = design_name)
}

trajectory_accuracy_all <- bind_rows(
  summarize_traj_accuracy(traj_targeted_nested, "Targeted", "Nested"),
  summarize_traj_accuracy(traj_targeted_fixed,  "Targeted", "Fixed Panel"),
  summarize_traj_accuracy(traj_protein_nested,  "Protein",  "Nested"),
  summarize_traj_accuracy(traj_protein_fixed,   "Protein",  "Fixed Panel"),
  summarize_traj_accuracy(traj_combined_nested, "Combined", "Nested"),
  summarize_traj_accuracy(traj_combined_fixed,  "Combined", "Fixed Panel")
)

cat("\nPercent-predicted-fatigued trajectory summary (all panels/designs):\n")
print(trajectory_accuracy_all %>% mutate(pct_correct = round(pct_correct, 3)), n = Inf)

write_csv(trajectory_accuracy_all, here("Results", "Fatigue_PercentPredicted_Trajectories_AllPanels.csv"))
cat("Saved: Results/Fatigue_PercentPredicted_Trajectories_AllPanels.csv\n")

#-------- Figure 3: Predicted-Fatigued Classification Rate (Nested) --------#

day_labels_df <- tibble(x = c(1.5, 4, 6.5), label = c("Day 1", "Day 2", "Day 3"))

df_acc_all3 <- trajectory_accuracy_all %>%
  filter(Design == "Nested") %>%
  mutate(Panel = factor(Panel, levels = c("Targeted", "Protein", "Combined")))

p_acc_all3 <- ggplot(df_acc_all3, aes(x = Measure, y = pct_correct, group = Panel)) +

  geom_rect(data = phase_shading,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
            inherit.aes = FALSE, alpha = 0.1) +
  scale_fill_manual(name = "Phase",
                     values = c("Rested" = "steelblue", "Fatigued" = "firebrick", "Recovery" = "forestgreen"),
                     breaks = c("Rested", "Fatigued", "Recovery")) +

  new_scale_fill() +

  geom_col(aes(fill = Panel), position = position_dodge(width = 0.7), alpha = 0.85, width = 0.65) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.7), linewidth = 0.35) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_vline(xintercept = day_boundaries, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(data = day_labels_df, aes(x = x, y = 1.10, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 3.5) +

  scale_fill_manual(name = "Model",
                     values = c("Targeted" = "#E69F00", "Protein" = "#0072B2", "Combined" = "#009E73"),
                     labels = c("Targeted" = "Targeted Biomarker Panel", "Protein" = "Protein Panel",
                                "Combined" = "Combined Panel")) +

  scale_x_discrete(labels = x_labels) +
  scale_y_continuous(limits = c(0, 1.15), expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +

  labs(x = NULL, y = "Predicted Fatigued (%)", title = "") +
  theme_classic() +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold", hjust = 0.5, size = 10))

ggsave(here("Results", "Figure3_PercentPredictedFatigued_Trajectories.tiff"),
       plot = p_acc_all3, width = 10, height = 6, units = "in", dpi = 300)
cat("Saved: Results/Figure3_PercentPredictedFatigued_Trajectories.tiff\n")

#=============================================================
# FIGURE 2: Temporal Changes in Panel Biomarkers (IgA, UA, ANG, C1S)
#=============================================================
# Shows the four biomarkers retained in the final combined panel,
# using the same visual style (phase shading, boxplot + individual
# trajectories + group mean line, Tukey-adjusted significance brackets for
# the 5 established contrasts) as the full 6-marker analysis above. The
# full mixed-effects results for all 6 targeted markers are reported in
# text (Biomarker_Posthoc_Contrasts.csv, above); this figure narrows the
# visual focus to the final combined panel only.

cat("\n=== FIGURE 2: Temporal Changes in Panel Biomarkers (IgA, UA, ANG, C1S) ===\n")

df_panel_targeted <- df_long %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, log_IgA, log_UA) %>%
  pivot_longer(cols = c(log_IgA, log_UA), names_to = "Biomarker", values_to = "Value") %>%
  mutate(Biomarker = recode(Biomarker, "log_IgA" = "IgA", "log_UA" = "UA"))

df_panel_c1s <- df_imputed %>%
  filter(Accession == "NP_958850.1") %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, Value = log_Abundance) %>%
  mutate(Biomarker = "C1S")

df_panel_ang <- df_imputed %>%
  filter(Accession == "NP_001091046.1") %>%
  mutate(ID = as.character(ID), Measure = as.character(Measure)) %>%
  dplyr::select(ID, Measure, Value = log_Abundance) %>%
  mutate(Biomarker = "ANG")

df_panel_biomarkers <- bind_rows(df_panel_targeted, df_panel_c1s, df_panel_ang) %>%
  mutate(
    Biomarker = factor(Biomarker, levels = c("IgA", "UA", "ANG", "C1S")),
    Measure   = factor(Measure, levels = c("1","2","3","4","5","6","7"))
  )

panel_biomarker_models <- lapply(levels(df_panel_biomarkers$Biomarker), function(b) {
  dat <- df_panel_biomarkers %>% filter(Biomarker == b) %>% mutate(Measure = relevel(Measure, ref = "3"))
  lmer(Value ~ Measure + (1 | ID), data = dat, REML = FALSE)
})
names(panel_biomarker_models) <- levels(df_panel_biomarkers$Biomarker)

panel_biomarker_posthoc <- lapply(names(panel_biomarker_models), function(b) {
  emm   <- emmeans(panel_biomarker_models[[b]], ~ Measure)
  pairs <- as.data.frame(pairs(emm, adjust = "tukey"))

  pairs %>%
    filter(contrast %in% c(
      "Measure1 - Measure3", "Measure3 - Measure1", "Measure1 - Measure2", "Measure2 - Measure1",
      "Measure3 - Measure4", "Measure4 - Measure3", "Measure4 - Measure5", "Measure5 - Measure4",
      "Measure4 - Measure6", "Measure6 - Measure4", "Measure4 - Measure7", "Measure7 - Measure4"
    )) %>%
    mutate(
      Biomarker = b,
      Pattern = case_when(
        contrast %in% c("Measure1 - Measure3", "Measure3 - Measure1") ~ "Baseline Stability",
        contrast %in% c("Measure1 - Measure2", "Measure2 - Measure1") ~ "Diurnal Effect",
        contrast %in% c("Measure3 - Measure4", "Measure4 - Measure3") ~ "Fatigue Response",
        contrast %in% c("Measure4 - Measure5", "Measure5 - Measure4") ~ "Immediate Recovery",
        TRUE ~ "Next Day Recovery"
      )
    )
}) %>% bind_rows()

panel_biomarker_posthoc_sig <- panel_biomarker_posthoc %>% filter(p.value < 0.05)

cat("\nSignificant post-hoc contrasts for panel biomarkers (IgA, UA, ANG, C1S):\n")
print(panel_biomarker_posthoc_sig %>%
        dplyr::select(Biomarker, Pattern, contrast, estimate, p.value) %>%
        mutate(across(where(is.numeric), ~round(.x, 3))))

write_csv(panel_biomarker_posthoc_sig %>%
            dplyr::select(Biomarker, Pattern, contrast, estimate, SE, p.value),
          here("Results", "Panel_Biomarker_Posthoc_Contrasts.csv"))

df_panel_biomarker_means <- df_panel_biomarkers %>%
  group_by(Biomarker, Measure) %>%
  summarise(mean_val = mean(Value, na.rm = TRUE), .groups = "drop")

pairs_plot_panel_biomarkers <- panel_biomarker_posthoc_sig %>%
  mutate(
    group1 = str_extract(contrast, "(?<=Measure)\\d(?= -)"),
    group2 = str_extract(contrast, "(?<=- Measure)\\d"),
    sig_label = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*"),
    Biomarker = factor(Biomarker, levels = c("IgA", "UA", "ANG", "C1S"))
  )

y_positions_panel_biomarkers <- df_panel_biomarkers %>%
  group_by(Biomarker) %>%
  summarise(y_max = max(Value, na.rm = TRUE), .groups = "drop")

pairs_plot_panel_biomarkers <- pairs_plot_panel_biomarkers %>%
  left_join(y_positions_panel_biomarkers, by = "Biomarker") %>%
  group_by(Biomarker) %>%
  mutate(y_position = y_max + (row_number() * 0.3)) %>%
  ungroup()

day_label_y <- pairs_plot_panel_biomarkers %>%
  group_by(Biomarker) %>%
  summarise(label_y = max(y_position, na.rm = TRUE) + 0.5, .groups = "drop")

day_labels_df_fig2 <- tidyr::crossing(
  Biomarker = levels(df_panel_biomarkers$Biomarker),
  tibble(x = c(1.5, 4, 6.5), label = c("Day 1", "Day 2", "Day 3"))
) %>%
  mutate(Biomarker = factor(Biomarker, levels = levels(df_panel_biomarkers$Biomarker))) %>%
  left_join(day_label_y, by = "Biomarker")

p_panel_biomarkers <- ggplot(df_panel_biomarkers, aes(x = Measure, y = Value)) +

  geom_rect(data = phase_shading,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = phase),
            inherit.aes = FALSE, alpha = 0.1) +
  scale_fill_manual(name = "Phase",
                     values = c("Rested" = "steelblue", "Fatigued" = "firebrick", "Recovery" = "forestgreen"),
                     breaks = c("Rested", "Fatigued", "Recovery")) +

  geom_boxplot(outlier.shape = NA, fill = "grey90", alpha = 0.5) +
  geom_line(aes(group = ID, color = ID), alpha = 0.4, linewidth = 0.5) +
  geom_point(aes(color = ID), alpha = 0.4, size = 1.5) +
  geom_line(data = df_panel_biomarker_means, aes(x = Measure, y = mean_val, group = 1),
            color = "black", linewidth = 1.2) +
  geom_point(data = df_panel_biomarker_means, aes(x = Measure, y = mean_val), color = "black", size = 3) +

  geom_signif(data = pairs_plot_panel_biomarkers,
              aes(xmin = group1, xmax = group2, annotations = sig_label, y_position = y_position),
              manual = TRUE, inherit.aes = FALSE, tip_length = 0.01, textsize = 3, vjust = 0.5) +

  geom_vline(xintercept = day_boundaries, linetype = "dashed", color = "grey40", linewidth = 0.5) +
  geom_text(data = day_labels_df_fig2, aes(x = x, y = label_y, label = label),
            inherit.aes = FALSE, fontface = "bold", size = 3.5) +

  facet_wrap(~ Biomarker, ncol = 2, scales = "free_y", axes = "all",
             labeller = as_labeller(c(
               "IgA" = str_wrap("Immunoglobulin A (IgA)", width = 30),
               "UA"  = str_wrap("Uric Acid (UA)", width = 30),
               "ANG" = str_wrap("Angiogenin (ANG)", width = 30),
               "C1S" = str_wrap("Complement C1s subcomponent isoform 1 preproprotein (C1S)", width = 30)
             ))) +

  scale_x_discrete(labels = x_labels) +

  labs(x = NULL, y = "Log Concentration / Abundance", color = "Participant", title = "") +
  theme_classic() +
  theme(legend.position = "bottom", strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9)) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))

ggsave(here("Results", "Figure2_Panel_Biomarker_Trajectories.tiff"),
       plot = p_panel_biomarkers, width = 11, height = 9, units = "in", dpi = 300)
cat("Saved: Results/Figure2_Panel_Biomarker_Trajectories.tiff\n")

#=============================================================
# END OF SCRIPT
#=============================================================

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
cat("\nFinal panels:\n")
cat("  Targeted:", paste(final_targeted_panel_names, collapse = " + "), "\n")
cat("  Protein: ", paste(final_protein_panel, collapse = " + "), "\n")
cat("  Combined:", paste(final_combined_panel, collapse = " + "), "\n")
cat("\nSee Results/ for all exported tables and figures.\n")
cat("Run PROSPER_Imputation_Sensitivity.R separately to reproduce the\n")
cat("imputation sensitivity analysis (Methods).\n")

plan(sequential)
