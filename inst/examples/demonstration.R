# ============================================================
# BayCAST — Demonstration
# ============================================================
# Self-contained semi-synthetic example.
#
#   - 400 subjects, 100 features
#   - 2 latent subgroups (cluster intercept gap = 0.4 on log time)
#   - 20 of 100 features are signal (rows 1..20 of the loading matrix)
#   - True beta = (-0.60, +0.40, +0.40) matching paper simulation
#   - Fits BayCAST in both modes:
#       PM  (sigma_spike2 = 0)    — point-mass spike-and-slab
#       CO  (sigma_spike2 = 0.01) — continuous narrow-normal spike
#   - Compares: K_+ recovery, variable-selection AUC, cluster ARI, beta.
#
# Runtime: ~10 seconds total on a laptop (PM + CO).
# Usage:   source(system.file("examples", "demonstration.R", package = "BayCAST"))
# ============================================================

suppressPackageStartupMessages({
  library(BayCAST)
})

set.seed(2026)

# ------------------------------------------------------------
# 1. Generate semi-synthetic data
# ------------------------------------------------------------
N         <- 400      # subjects
P         <- 100      # features
d_true    <- 3        # latent dimension
K_true    <- 2        # latent subgroups
n_signal  <- 20       # signal features
sigma_obs <- 0.20     # residual SD on log time

# Latent factor structure with 2 clusters
pi_true <- c(0.60, 0.40)
z_true  <- sample(1:K_true, N, replace = TRUE, prob = pi_true)
mu_t    <- c(4.0, 4.4)                                 # cluster intercepts on log time (gap = 0.4)
nu_U    <- list(c(-1, 0, 0), c(+1, 0, 0))             # cluster centres in factor space
U       <- t(sapply(z_true, function(k)
                    rnorm(d_true, mean = nu_U[[k]], sd = 0.5)))

# Sparse loading matrix: first n_signal rows are signal
A <- matrix(0, P, d_true)
signal_idx <- 1:n_signal
A[signal_idx, ] <- matrix(rnorm(n_signal * d_true, sd = 0.5), n_signal, d_true)

# Observed feature matrix
X <- U %*% t(A) + matrix(rnorm(N * P, sd = 0.3), N, P)

# Clinical covariates and true regression effect
# Same beta as paper simulation: (HXMI, SEX, RACE_B)
W <- cbind(
  HXMI   = rbinom(N, 1, 0.35),
  SEX    = rbinom(N, 1, 0.50),
  RACE_B = rbinom(N, 1, 0.25)
)
beta_true <- c(HXMI = -0.20, SEX = +0.10, RACE_B = +0.10)

# Latent log event time and current-status indicator
log_T  <- as.numeric(W %*% beta_true + mu_t[z_true] + rnorm(N, sd = sigma_obs))
log_C  <- log(runif(N, min = 25, max = 95))            # ages 25..95, log scale
status <- as.integer(log_T <= log_C)

cat(strrep("=", 60), "\n", sep = "")
cat(" BayCAST demo — semi-synthetic data\n")
cat(strrep("=", 60), "\n", sep = "")
cat(sprintf("  N = %d, P = %d, d_true = %d, K_true = %d\n",
            N, P, d_true, K_true))
cat(sprintf("  Signal features: rows 1..%d (%d of %d, %.0f%%)\n",
            n_signal, n_signal, P, 100 * n_signal / P))
cat(sprintf("  Event rate    : %.2f\n", mean(status)))
cat(sprintf("  True beta     : HXMI=%+.2f  SEX=%+.2f  RACE_B=%+.2f\n",
            beta_true[1], beta_true[2], beta_true[3]))
cat(sprintf("                  (matches paper simulation DGP)\n\n"))

# ------------------------------------------------------------
# 2. Fit PM
# ------------------------------------------------------------
cat(strrep("-", 60), "\n", sep = "")
cat(" [Fit 1] PM  (sigma_spike2 = 0)\n")
cat(strrep("-", 60), "\n", sep = "")
t0 <- Sys.time()
fit_pm <- baycast(
  time       = log_C,
  status     = status,
  X          = X,
  covariates = W,
  d          = d_true,
  iters      = 3000, burn = 1500, thin = 2,
  pi0        = 0.10,
  sigma_spike2 = 0.0,
  seed       = 1L,
  verbose    = FALSE
)
cat(sprintf("  elapsed : %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ------------------------------------------------------------
# 3. Fit CO
# ------------------------------------------------------------
cat("\n", strrep("-", 60), "\n", sep = "")
cat(" [Fit 2] CO  (sigma_spike2 = 0.01)\n")
cat(strrep("-", 60), "\n", sep = "")
t0 <- Sys.time()
fit_co <- baycast(
  time       = log_C,
  status     = status,
  X          = X,
  covariates = W,
  d          = d_true,
  iters      = 3000, burn = 1500, thin = 2,
  pi0        = 0.10,
  sigma_spike2 = 0.01,
  seed       = 1L,
  verbose    = FALSE
)
cat(sprintf("  elapsed : %.1f s\n",
            as.numeric(difftime(Sys.time(), t0, units = "secs"))))

# ------------------------------------------------------------
# 4. Recovery summary
# ------------------------------------------------------------
auc_score <- function(truth, score) {
  # Higher score should be more likely to be signal (truth = 1).
  # Mann-Whitney U formulation.
  r     <- rank(score)
  n_pos <- sum(truth == 1); n_neg <- sum(truth == 0)
  (sum(r[truth == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

modal <- function(x) as.integer(names(sort(table(x), decreasing = TRUE))[1])

summarise_fit <- function(fit) {
  pip <- fit$pip
  zd  <- fit$z_draws
  z_mode <- if (nrow(zd) == N) apply(zd, 1, modal) else apply(zd, 2, modal)

  truth_vec <- rep(0L, P); truth_vec[signal_idx] <- 1L
  list(
    K_mode      = modal(fit$Kplus_draws),
    PIP_gt_0.5  = sum(pip > 0.5),
    varsel_AUC  = auc_score(truth_vec, pip),
    cluster_ARI = if (requireNamespace("mclust", quietly = TRUE))
                    mclust::adjustedRandIndex(z_mode, z_true)
                  else NA_real_,
    beta_mean   = rowMeans(fit$beta_draws),
    beta_lo     = apply(fit$beta_draws, 1, quantile, 0.025),
    beta_hi     = apply(fit$beta_draws, 1, quantile, 0.975)
  )
}

res_pm <- summarise_fit(fit_pm)
res_co <- summarise_fit(fit_co)

cat("\n", strrep("=", 60), "\n", sep = "")
cat(" RECOVERY COMPARISON: PM vs CO\n")
cat(strrep("=", 60), "\n", sep = "")
cat(sprintf("%-32s  %10s  %10s\n", "", "PM", "CO"))
cat(sprintf("%-32s  %10s  %10s\n",
            strrep("-", 32), strrep("-", 10), strrep("-", 10)))
cat(sprintf("%-32s  %10d  %10d   (truth %d)\n",
            "Posterior modal K_+",
            res_pm$K_mode, res_co$K_mode, K_true))
cat(sprintf("%-32s  %10d  %10d   (truth %d)\n",
            "# features with PIP > 0.5",
            res_pm$PIP_gt_0.5, res_co$PIP_gt_0.5, n_signal))
cat(sprintf("%-32s  %10.3f  %10.3f\n",
            "Variable-selection AUC",
            res_pm$varsel_AUC, res_co$varsel_AUC))
cat(sprintf("%-32s  %10.3f  %10.3f\n",
            "Cluster ARI (vs truth)",
            res_pm$cluster_ARI, res_co$cluster_ARI))
cat(sprintf("%-32s  %+10.3f  %+10.3f   (truth %+.2f)\n",
            "beta_HXMI   posterior mean",
            res_pm$beta_mean[1], res_co$beta_mean[1], beta_true[1]))
cat(sprintf("%-32s  %+10.3f  %+10.3f   (truth %+.2f)\n",
            "beta_SEX    posterior mean",
            res_pm$beta_mean[2], res_co$beta_mean[2], beta_true[2]))
cat(sprintf("%-32s  %+10.3f  %+10.3f   (truth %+.2f)\n",
            "beta_RACE_B posterior mean",
            res_pm$beta_mean[3], res_co$beta_mean[3], beta_true[3]))

# ------------------------------------------------------------
# 5. Visualisation — PIP and cluster posterior similarity
# ------------------------------------------------------------
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))

# 5a. PIP histograms (PM vs CO)
hist(fit_pm$pip, breaks = 40, col = "#7AA6D9", border = "white",
     xlim = c(0, 1), main = "PM: posterior inclusion probability",
     xlab = "PIP")
abline(v = 0.5, col = "red", lty = 2)
rug(fit_pm$pip[signal_idx], col = "darkgreen", lwd = 2)

hist(fit_co$pip, breaks = 40, col = "#9FD17A", border = "white",
     xlim = c(0, 1), main = "CO: posterior inclusion probability",
     xlab = "PIP")
abline(v = 0.5, col = "red", lty = 2)
rug(fit_co$pip[signal_idx], col = "darkgreen", lwd = 2)

# 5b. PIP scatter (signal in green, noise in grey)
col_pt <- ifelse(seq_len(P) %in% signal_idx, "darkgreen", "grey50")
plot(fit_pm$pip, fit_co$pip, pch = 19, col = col_pt,
     xlab = "PM PIP", ylab = "CO PIP",
     xlim = c(0, 1), ylim = c(0, 1),
     main = "PM PIP vs CO PIP (signal = green)")
abline(0, 1, lty = 2, col = "grey60")
abline(h = 0.5, v = 0.5, lty = 3, col = "red")

# 5c. Kplus posterior bar plots — use common range
kp_range <- range(c(fit_pm$Kplus_draws, fit_co$Kplus_draws))
levels_kp <- seq(kp_range[1], kp_range[2])
kp <- rbind(
  PM = as.numeric(table(factor(fit_pm$Kplus_draws, levels = levels_kp))),
  CO = as.numeric(table(factor(fit_co$Kplus_draws, levels = levels_kp)))
)
colnames(kp) <- levels_kp
barplot(kp, beside = TRUE, col = c("#7AA6D9", "#9FD17A"),
        main = "Posterior of K_+   (truth = 2)",
        xlab = "K_+", ylab = "Frequency",
        legend.text = c("PM", "CO"),
        args.legend = list(x = "topright", bty = "n"))

mtext("BayCAST demo  —  PM (tau²=0) vs CO (tau²=0.01)",
      outer = TRUE, cex = 1.1, font = 2)
par(op)

cat("\n[DONE]\n")