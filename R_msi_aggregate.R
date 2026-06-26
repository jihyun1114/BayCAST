# R_msi_aggregate.R
# Aggregate 18 task summaries into paper-style tables.
# Run after rsync from MSI.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

ROOT <- "output/msi_grid_full"

# ----------------------------------------------------------------
# 1. Load all
# ----------------------------------------------------------------
all_csv <- list.files(ROOT, pattern = "summary\\.csv$",
                      recursive = TRUE, full.names = TRUE)
cat("Found", length(all_csv), "summary files\n\n")

big <- do.call(rbind, lapply(all_csv, read.csv, stringsAsFactors = FALSE))
cat("Total rows:", nrow(big), " (expect 100 reps × 18 = 1800)\n")
cat("Per scenario × mode:\n")
print(table(big$scenario, big$mode))

write.csv(big, file.path(ROOT, "ALL_REPS.csv"), row.names = FALSE)
cat("\n[Saved] ALL_REPS.csv\n\n")

# ----------------------------------------------------------------
# 2. Paper Table 2 — Variable selection AUC
# ----------------------------------------------------------------
cat(strrep("=", 70), "\n")
cat(" Table 2 — Variable selection AUC (proposed)\n")
cat(strrep("=", 70), "\n")
tab_varsel <- aggregate(varsel_auc ~ scenario + mode,
                        data = big,
                        FUN = function(x) c(mean = mean(x), sd = sd(x)))
tab_varsel <- do.call(data.frame, tab_varsel)
colnames(tab_varsel) <- c("scenario", "mode", "mean", "sd")
tab_varsel$value <- sprintf("%.3f (%.3f)", tab_varsel$mean, tab_varsel$sd)
print(tab_varsel[, c("scenario", "mode", "value")], row.names = FALSE)
write.csv(tab_varsel, file.path(ROOT, "Table2_varsel_AUC.csv"), row.names = FALSE)

# ----------------------------------------------------------------
# 3. Paper Table 3 — Mean K+ estimate
# ----------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat(" Table 3 — Mean K+ posterior mode (proposed)\n")
cat(strrep("=", 70), "\n")
tab_kplus <- aggregate(cbind(Kp_mode, true_K) ~ scenario + mode,
                       data = big, FUN = mean)
tab_kplus$Kp_mode <- round(tab_kplus$Kp_mode, 2)
print(tab_kplus, row.names = FALSE)
write.csv(tab_kplus, file.path(ROOT, "Table3_kplus.csv"), row.names = FALSE)

# ----------------------------------------------------------------
# 4. Paper Table 4 — Test ARI
# ----------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat(" Table 4 — Test ARI (proposed)\n")
cat(strrep("=", 70), "\n")
tab_ari <- aggregate(test_ari ~ scenario + mode,
                     data = big,
                     FUN = function(x) c(mean = mean(x, na.rm=TRUE),
                                          sd = sd(x, na.rm=TRUE)))
tab_ari <- do.call(data.frame, tab_ari)
colnames(tab_ari) <- c("scenario", "mode", "mean", "sd")
tab_ari$value <- sprintf("%.3f (%.3f)", tab_ari$mean, tab_ari$sd)
print(tab_ari[, c("scenario", "mode", "value")], row.names = FALSE)
write.csv(tab_ari, file.path(ROOT, "Table4_test_ARI.csv"), row.names = FALSE)

# ----------------------------------------------------------------
# 5. Paper Table 5 — Test Y-AUC
# ----------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat(" Table 5 — Test Y-AUC (proposed)\n")
cat(strrep("=", 70), "\n")
tab_yauc <- aggregate(test_yauc ~ scenario + mode,
                      data = big,
                      FUN = function(x) c(mean = mean(x, na.rm=TRUE),
                                           sd = sd(x, na.rm=TRUE)))
tab_yauc <- do.call(data.frame, tab_yauc)
colnames(tab_yauc) <- c("scenario", "mode", "mean", "sd")
tab_yauc$value <- sprintf("%.3f (%.3f)", tab_yauc$mean, tab_yauc$sd)
print(tab_yauc[, c("scenario", "mode", "value")], row.names = FALSE)
write.csv(tab_yauc, file.path(ROOT, "Table5_test_YAUC.csv"), row.names = FALSE)

# ----------------------------------------------------------------
# 6. Paper Table 6 — β inference (Bias, ASD, ESE, CP)
# ----------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat(" Table 6 — β inference (proposed)\n")
cat(strrep("=", 70), "\n")

COV_NAMES <- c("HXMI", "SEX", "RACE_B")
beta_rows <- list()

for (s in unique(big$scenario)) {
  for (m in unique(big$mode)) {
    sub <- big[big$scenario == s & big$mode == m, ]
    if (nrow(sub) == 0) next
    row <- c(scenario = s, mode = m)
    for (cv in COV_NAMES) {
      bias <- mean(sub[[paste0("beta_bias_", cv)]])
      asd  <- mean(sub[[paste0("beta_sd_",   cv)]])
      ese  <- sd  (sub[[paste0("beta_hat_",  cv)]])
      cp   <- mean(sub[[paste0("beta_cov_",  cv)]])
      row <- c(row,
               setNames(c(round(bias, 3), round(asd, 3),
                          round(ese, 3), round(cp, 2)),
                        paste0(cv, c("_bias", "_ASD", "_ESE", "_CP"))))
    }
    beta_rows[[length(beta_rows) + 1]] <- row
  }
}
tab_beta <- as.data.frame(do.call(rbind, beta_rows), stringsAsFactors = FALSE)
print(tab_beta, row.names = FALSE)
write.csv(tab_beta, file.path(ROOT, "Table6_beta_inference.csv"), row.names = FALSE)

cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("\n[Saved] all paper-style tables to %s/\n", ROOT))
cat("  - Table2_varsel_AUC.csv\n")
cat("  - Table3_kplus.csv\n")
cat("  - Table4_test_ARI.csv\n")
cat("  - Table5_test_YAUC.csv\n")
cat("  - Table6_beta_inference.csv\n")
cat("  - ALL_REPS.csv (raw rep-level)\n")
cat("[DONE]\n")
