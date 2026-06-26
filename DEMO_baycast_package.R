# ============================================================
# DEMO_baycast_package.R
#
# 미팅 시연용 — baycast 패키지가 CATHGEN 실데이터로 실제로 돌아간다
# 는 것을 한 화면에 보여주는 스크립트.
#
# Usage (run from project root ~/baycast):
#   source("DEMO_baycast_package.R")
#
# Two modes:
#   QUICK = TRUE   : iters=500,  burn=250,  ~1.5 min   ← 미팅 중 라이브용
#   QUICK = FALSE  : iters=4000, burn=2000, ~10 min    ← 미리 돌려놓는 용
# ============================================================

QUICK <- TRUE

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" baycast package — live demo on CATHGEN data\n")
cat(strrep("=", 60), "\n\n", sep = "")

# ============================================================
# 1. Load the installed package
# ============================================================
library(BayCAST)
cat("Package version:", as.character(packageVersion("baycast")), "\n")
cat("Source: github.com/jihyun1114/baycast\n\n")

# ============================================================
# 2. Load CATHGEN preprocessed data
#    (already cleaned in earlier R01 step)
# ============================================================
stopifnot(dir.exists("output/01_preprocessing"))

X    <- as.matrix(read.csv("output/01_preprocessing/GeneCovariates.csv",
                           row.names = 1, check.names = FALSE))
clin <- read.csv("output/01_preprocessing/Covariates.csv",
                 row.names = 1, check.names = FALSE)
clin <- clin[rownames(X), , drop = FALSE]

# survival outcome + forced covariates
time   <- as.numeric(clin$AGE)
status <- as.integer(clin$cadIndex_binary)

W <- cbind(
  HXMI   = as.integer(clin$HXMI),
  SEX    = as.integer(clin$SEX),
  RACE_B = as.integer(clin$RACE == 2)
)
W[is.na(W)] <- 0L
rownames(W) <- rownames(X)

cat("Data loaded:\n")
cat(sprintf("  N subjects : %d\n", nrow(X)))
cat(sprintf("  P features : %d\n", ncol(X)))
cat(sprintf("  Event rate : %.3f\n", mean(status)))
cat(sprintf("  Covariates : %s\n", paste(colnames(W), collapse = ", ")))
cat("\n")

# ============================================================
# 3. Fit
# ============================================================
if (QUICK) {
  iters <- 500L; burn <- 250L; thin <- 5L
  cat("MODE: QUICK demo (iters=500). For paper-quality fit set QUICK = FALSE.\n\n")
} else {
  iters <- 4000L; burn <- 2000L; thin <- 10L
  cat("MODE: FULL fit (iters=4000) — matches main run config.\n\n")
}

t0 <- Sys.time()

fit <- baycast(
  time       = time,
  status     = status,
  X          = X,
  covariates = W,
  d          = 9,
  K_init     = 15,
  K_max      = 50,
  iters      = iters,
  burn       = burn,
  thin       = thin,
  pi0        = 0.02,
  sigmaA2    = 1.0,
  slice      = "kgw",
  init_A     = "pca",
  init_z     = "kmeans",
  seed       = 2026,
  verbose    = TRUE
)

elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
cat(sprintf("\n[FIT DONE] elapsed = %.2f min\n\n", elapsed))

# ============================================================
# 4. Show results
# ============================================================
cat(strrep("-", 60), "\n", sep = "")
cat(" print(fit)\n")
cat(strrep("-", 60), "\n", sep = "")
print(fit)

cat("\n", strrep("-", 60), "\n", sep = "")
cat(" summary(fit)\n")
cat(strrep("-", 60), "\n\n", sep = "")
print(summary(fit))

# ============================================================
# 5. Cluster-level breakdown
# ============================================================
cat("\n", strrep("-", 60), "\n", sep = "")
cat(" Cluster-level event rate (posterior-mode partition)\n")
cat(strrep("-", 60), "\n", sep = "")

z_mode <- apply(fit$z_draws, 1, function(zi) {
  tab <- table(zi); as.integer(names(which.max(tab)))
})
br <- table(z_mode)
for (k in names(br)) {
  in_k <- z_mode == as.integer(k)
  cat(sprintf("  cluster %s : n = %3d   event rate = %.3f\n",
              k, sum(in_k), mean(status[in_k])))
}

# ============================================================
# 6. Top selected genes
# ============================================================
cat("\n", strrep("-", 60), "\n", sep = "")
cat(" Top 10 features by posterior inclusion probability\n")
cat(strrep("-", 60), "\n", sep = "")

loading_norm <- sqrt(rowSums(fit$A_mean^2))
top_idx <- order(fit$pip, loading_norm, decreasing = TRUE)[1:10]
top_tab <- data.frame(
  rank        = 1:10,
  feature     = colnames(X)[top_idx],
  PIP         = round(fit$pip[top_idx], 4),
  loading_norm = round(loading_norm[top_idx], 3)
)
print(top_tab, row.names = FALSE)

# ============================================================
# 7. Headline numbers
# ============================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat(" HEADLINE\n")
cat(strrep("=", 60), "\n", sep = "")

beta_mean <- rowMeans(fit$beta_draws)
beta_ci   <- t(apply(fit$beta_draws, 1, quantile, c(0.025, 0.975)))

cat(sprintf("  K+ posterior mode               : %d\n",
            as.integer(names(sort(table(fit$Kplus_draws), decreasing = TRUE))[1])))
cat(sprintf("  Cluster sizes (mode partition)  : %s\n",
            paste(sort(as.integer(table(z_mode)), decreasing = TRUE), collapse = " / ")))
cat(sprintf("  PIP > 0.5  selected features    : %d / %d\n",
            sum(fit$pip > 0.5), length(fit$pip)))
cat(sprintf("  PIP > 0.95 selected features    : %d / %d\n",
            sum(fit$pip > 0.95), length(fit$pip)))
for (j in seq_along(beta_mean)) {
  cat(sprintf("  beta_%-8s  = %+.3f  95%%CI [%+.3f, %+.3f]\n",
              colnames(W)[j], beta_mean[j],
              beta_ci[j, 1], beta_ci[j, 2]))
}

# ============================================================
# 8. Save fit object
# ============================================================
out_path <- if (QUICK) "demo_fit_quick.rds" else "demo_fit_full.rds"
saveRDS(fit, out_path)
cat(sprintf("\n[Saved] %s\n", out_path))

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Demo complete.\n")
cat(strrep("=", 60), "\n", sep = "")