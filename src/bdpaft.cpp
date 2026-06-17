// bdpaft_continuous_ssvs.cpp
// BDPAFT sampler with CONTINUOUS spike-slab variable selection.
//
// Difference from point-mass version (bdpaft.cpp):
//   point-mass:  A_{j.} | delta_j = 0  has mass at exactly 0
//                A_{j.} | delta_j = 1  ~ N(0, sigmaA2 * I_d)   (slab)
//
//   continuous:  A_{j.} | delta_j = 0  ~ N(0, tau_spike2 * I_d)  (narrow spike)
//                A_{j.} | delta_j = 1  ~ N(0, sigmaA2     * I_d)  (slab)
//
// Activated by passing tau_spike2 > 0.
// Setting tau_spike2 <= 0 reproduces the original point-mass behavior.
//
// Recommended values:
//   tau_spike2 / sigmaA2  in  [1e-4, 1e-2]
//   e.g. sigmaA2 = 1.0, tau_spike2 = 0.01 (spike SD = 0.1 vs slab SD = 1.0)
//
// MCMC: Gibbs sampling with Kalli-Griffin-Walker slice for DP truncation.
// Model:
//   log T_i = W_i'beta + mu_{z_i} + epsilon_i,  epsilon ~ N(0, sigma^2_{z_i})
//   z_i ~ DP(alpha, G_0)         (cluster assignment, KGW slice)
//   X_i = U_i A' + e_i           (factor model with continuous spike-slab on A)
//   U_i | z_i ~ N(nu_{z_i}, Sigma_{z_i})
//
// [[Rcpp::plugins(cpp17)]]
// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iomanip>

using namespace Rcpp;
using namespace arma;

// ========================== utilities ==========================

inline double clamp01(double x, double lo, double hi){
  if(x < lo) return lo;
  if(x > hi) return hi;
  return x;
}

inline double dnorm_log(double x, double m, double s){
  const double c = -0.5 * std::log(2.0 * arma::datum::pi);
  double z = (x - m) / s;
  return c - std::log(s) - 0.5 * z * z;
}

inline arma::vec rmvnorm(const arma::vec& mu, const arma::mat& Sigma){
  arma::mat L = arma::chol(Sigma, "lower");
  return mu + L * arma::randn<arma::vec>(mu.n_elem);
}

inline arma::mat inv_sympd_ridge(const arma::mat& M,
                                 double ridge0 = 1e-10, int max_try = 8){
  arma::mat out; double ridge = ridge0;
  for(int t = 0; t < max_try; ++t){
    arma::mat A = M + ridge * arma::eye<arma::mat>(M.n_rows, M.n_cols);
    if(arma::inv_sympd(out, A)) return out;
    ridge *= 10.0;
  }
  return arma::pinv(M);
}

inline double logdet_ridge(const arma::mat& M,
                           double ridge0 = 1e-10, int max_try = 8){
  double sign = 0.0, ld = 0.0, ridge = ridge0;
  for(int t = 0; t < max_try; ++t){
    arma::mat A = M + ridge * arma::eye<arma::mat>(M.n_rows, M.n_cols);
    if(arma::log_det(ld, sign, A) && std::isfinite(ld)) return ld;
    ridge *= 10.0;
  }
  arma::vec ev = arma::clamp(arma::eig_sym(
    M + 1e-6 * arma::eye<arma::mat>(M.n_rows, M.n_cols)), 1e-15, arma::datum::inf);
  return arma::accu(arma::log(ev));
}

inline double rtruncnorm_right_icdf(double m, double s, double b){
  double Phi = clamp01(R::pnorm((b-m)/s, 0.0, 1.0, 1, 0), 1e-15, 1.0-1e-15);
  double u   = clamp01(R::runif(0.0, Phi), 1e-15, 1.0-1e-15);
  return m + s * R::qnorm(u, 0.0, 1.0, 1, 0);
}

inline double rtruncnorm_left_icdf(double m, double s, double a){
  double Phi = clamp01(R::pnorm((a-m)/s, 0.0, 1.0, 1, 0), 1e-15, 1.0-1e-15);
  double u   = clamp01(R::runif(Phi, 1.0), 1e-15, 1.0-1e-15);
  return m + s * R::qnorm(u, 0.0, 1.0, 1, 0);
}

inline arma::mat rinvwishart(const arma::mat& Psi, int df){
  return inv_sympd_ridge(arma::wishrnd(inv_sympd_ridge(Psi), df));
}

// ========================== init helpers ==========================

static void pca_init_X(const arma::mat& X, int d,
                       arma::mat& U_init, arma::mat& A_init, arma::vec& psi_init){
  int n = (int)X.n_rows, p = (int)X.n_cols;
  arma::mat Xc = X.each_row() - arma::mean(X, 0);
  arma::mat U, V; arma::vec s;
  arma::svd_econ(U, s, V, Xc, "both");
  int d_use = std::min(d, (int)s.n_elem);
  U_init = U.cols(0, d_use-1) * arma::diagmat(s.subvec(0, d_use-1));
  A_init.zeros(p, d);
  A_init.cols(0, d_use-1) = V.cols(0, d_use-1);
  arma::mat Rm = Xc - U_init * A_init.t();
  psi_init.set_size(p);
  for(int j = 0; j < p; ++j){
    double v = arma::as_scalar(arma::mean(arma::square(Rm.col(j))));
    psi_init(j) = (std::isfinite(v) && v > 1e-8) ? v : 1e-8;
  }
  if(d_use < d){
    arma::mat pad(n, d, fill::zeros);
    pad.cols(0, d_use-1) = U_init;
    for(int k = d_use; k < d; ++k) pad.col(k) = 0.1 * arma::randn<arma::vec>(n);
    U_init = pad;
  }
}

static void init_z_kmeans_helper(const arma::mat& U_pca, int K_init, arma::ivec& z){
  int n = (int)U_pca.n_rows;
  int K_use = std::min(K_init, n);
  arma::mat centroids;
  bool ok = arma::kmeans(centroids, U_pca.t(),
                         (arma::uword)K_use, arma::random_subset, 20, false);
  z.set_size(n);
  if(!ok || centroids.has_nan() || centroids.has_inf()){
    arma::uvec idx = arma::randperm((arma::uword)n);
    for(int i = 0; i < n; ++i) z((int)idx((arma::uword)i)) = i % K_init;
    return;
  }
  for(int i = 0; i < n; ++i){
    arma::vec xi = U_pca.row(i).t();
    arma::vec dists(K_use);
    for(int k = 0; k < K_use; ++k){
      arma::vec diff = xi - centroids.col(k);
      dists(k) = arma::dot(diff, diff);
    }
    z(i) = (int)dists.index_min();
  }
}

static void init_z_balanced_random_helper(int n, int K_init, arma::ivec& z){
  arma::uvec idx = arma::randperm((arma::uword)n);
  z.set_size(n);
  for(int i = 0; i < n; ++i) z((int)idx((arma::uword)i)) = i % K_init;
}

// ========================== alpha update (Escobar-West) ==========================

static void update_alpha_ew(double& alpha, int n, int K_plus,
                             double a_alpha, double b_alpha){
  double eta     = R::rbeta(alpha + 1.0, (double)n);
  double log_eta = std::log(std::max(eta, 1e-15));
  double denom   = std::max((double)n*(b_alpha-log_eta)+a_alpha+K_plus-1.0, 1e-12);
  double w       = (a_alpha+K_plus-1.0) / denom;
  double shape   = (R::runif(0.0,1.0) < w) ? (a_alpha+K_plus) : (a_alpha+K_plus-1.0);
  double scale   = 1.0 / std::max(b_alpha-log_eta, 1e-12);
  alpha = R::rgamma(std::max(shape, 1e-6), scale);
  if(!std::isfinite(alpha) || alpha <= 1e-6) alpha = 1e-6;
}

static void update_a_hyperprior(double& a_alpha, double alpha_dp, double b_alpha,
                                  double alpha_a_prior, double beta_a_prior,
                                  double proposal_sd,
                                  int& n_accept, int& n_total){
  if(alpha_dp <= 0.0 || b_alpha <= 0.0) return;
  double log_a_curr = std::log(std::max(a_alpha, 1e-12));
  double log_a_prop = log_a_curr + R::rnorm(0.0, proposal_sd);
  double a_prop     = std::exp(log_a_prop);

  auto logp = [&](double aa, double log_aa) {
    return aa * std::log(b_alpha) - std::lgamma(aa) + (aa - 1.0) * std::log(alpha_dp)
         + (alpha_a_prior - 1.0) * std::log(aa) - beta_a_prior * aa
         + log_aa;
  };
  double lp_curr = logp(a_alpha, log_a_curr);
  double lp_prop = logp(a_prop,  log_a_prop);

  n_total++;
  double log_accept = lp_prop - lp_curr;
  if(std::isfinite(log_accept) && std::log(R::runif(0.0,1.0)) < log_accept){
    a_alpha = a_prop;
    n_accept++;
  }
  if(!std::isfinite(a_alpha) || a_alpha <= 1e-6) a_alpha = 1e-6;
}

static void update_b_hyperprior(double& b_alpha, double alpha_dp, double a_alpha,
                                  double alpha_b_prior, double beta_b_prior){
  if(alpha_dp <= 0.0 || a_alpha <= 0.0) return;
  double shape = a_alpha + alpha_b_prior;
  double rate  = alpha_dp + beta_b_prior;
  b_alpha = R::rgamma(std::max(shape, 1e-6), 1.0/std::max(rate, 1e-12));
  if(!std::isfinite(b_alpha) || b_alpha <= 1e-6) b_alpha = 1e-6;
}

// ============================================================================
// [[Rcpp::export]]
Rcpp::List bdpaft_cpp(
  const arma::vec&     Y,
  const arma::vec&     logC,
  const arma::mat&     X,
  const arma::mat&     W,
  const arma::vec&     b0_beta,
  const arma::mat&     B0_beta,
  double               mu0,
  double               kappa_mu0,
  double               a0,
  double               b0,
  int                  d,
  int                  K_max,
  int                  K_init,
  int                  iters,
  int                  burn,
  int                  thin,
  double               sigmaA2,
  double               tau_spike2,   // NEW: <=0 -> point-mass spike, >0 -> continuous spike
  double               pi0,
  double               alpha_fixed,
  double               alpha_init,
  double               a_alpha,
  double               b_alpha,
  bool                 sample_ab,
  double               alpha_a_prior,
  double               beta_a_prior,
  double               alpha_b_prior,
  double               beta_b_prior,
  double               mh_a_sd,
  bool                 init_A_pca,
  bool                 init_z_kmeans,
  bool                 ind_slice,
  double               rho,
  double               nu_w_prior,
  const IntegerVector& diag_feat_idx,
  int                  diag_max_keep
){
  const int n = X.n_rows, p = X.n_cols;
  const int q = W.n_cols;
  if((int)Y.n_elem != n || (int)logC.n_elem != n) stop("Y,logC length must equal nrow(X)");
  if(q > 0 && (int)W.n_rows != n)                stop("W must have same nrow as X");
  if(q > 0 && ((int)b0_beta.n_elem != q ||
               (int)B0_beta.n_rows != q ||
               (int)B0_beta.n_cols != q))        stop("b0_beta/B0_beta dims must match ncol(W)");
  if(d < 1)                          stop("d must be >= 1");
  if(K_max < 2)                      stop("K_max must be >= 2");
  if(K_init < 2 || K_init > K_max)   stop("K_init must be in [2, K_max]");
  if(iters <= burn)                  stop("iters must be > burn");
  if(thin <= 0)                      stop("thin must be > 0");
  if(diag_max_keep < 1) diag_max_keep = 1;

  if(init_z_kmeans && !init_A_pca){
    Rcpp::warning("init_z_kmeans=TRUE requires init_A_pca=TRUE; falling back to balanced_random for z");
    init_z_kmeans = false;
  }

  pi0     = clamp01(pi0, 1e-8, 1.0-1e-8);
  sigmaA2 = std::max(sigmaA2, 1e-12);

  // NEW: continuous-SSVS toggle
  //   tau_spike2 <= 0     -> point-mass spike (legacy behavior)
  //   0 < tau_spike2 < .. -> narrow Normal spike (continuous SSVS)
  const bool use_continuous_ssvs = (tau_spike2 > 0.0);
  if(use_continuous_ssvs){
    if(tau_spike2 >= sigmaA2){
      Rcpp::warning("tau_spike2 should be < sigmaA2; clamping to sigmaA2/100");
      tau_spike2 = sigmaA2 / 100.0;
    }
    tau_spike2 = std::max(tau_spike2, 1e-12);
  }

  if(ind_slice) rho = clamp01(rho, 0.0, 1.0-1e-4);

  const bool sample_alpha = (alpha_fixed <= 0.0);
  double alpha = sample_alpha ? std::max(alpha_init, 1e-6) : alpha_fixed;
  if(!std::isfinite(alpha) || alpha <= 0) alpha = 1.0;

  int a_accept_ct = 0, a_total_ct = 0;

  const double logit_pi0 = std::log(pi0 / (1.0-pi0));
  const double t_eps = 0.5;

  // ---- diagnostic indices ----
  std::vector<int> diag0;
  for(int i = 0; i < diag_feat_idx.size(); ++i){
    int jj = diag_feat_idx[i];
    if(jj >= 1 && jj <= p) diag0.push_back(jj-1);
  }
  std::sort(diag0.begin(), diag0.end());
  diag0.erase(std::unique(diag0.begin(), diag0.end()), diag0.end());
  const int J = (int)diag0.size();
  if(J < 1) stop("diag_feat_idx has no valid indices in 1..p.");

  // ========================== Hyperparameters ==========================
  const double apsi = 2.0, bpsi = 2.0;
  const arma::vec nu0  = arma::zeros(d);
  const double kappa0  = 1.0;
  double nu_w = nu_w_prior;
  if(nu_w <= (double)d + 1.0){
    Rcpp::warning("nu_w_prior must be > d+1; clamping to d+2");
    nu_w = (double)d + 2.0;
  }
  double Psi0_scale = -1.0;
  bool   auto_psi0  = true;

  arma::mat B0_inv;
  arma::vec B0inv_b0;
  if(q > 0){
    B0_inv   = inv_sympd_ridge(B0_beta);
    B0inv_b0 = B0_inv * b0_beta;
  }

  // ========================== State ==========================
  arma::vec  t(n, fill::zeros);
  arma::ivec z(n, fill::zeros);
  arma::vec  s(n, fill::zeros);
  arma::mat  U(n, d, fill::randn);
  arma::mat  A(p, d, fill::zeros);
  arma::vec  psi(p, fill::ones);
  arma::vec  delta(p, fill::ones);

  arma::vec  beta_cur = (q > 0) ? b0_beta : arma::vec();

  arma::vec  V(K_max, fill::zeros);
  arma::vec  w_vec(K_max, fill::zeros);
  arma::vec  mu(K_max, fill::zeros);
  arma::vec  sig2(K_max, fill::ones);
  arma::mat  nu_mat(K_max, d, fill::zeros);
  arma::cube Sigma(d, d, K_max, fill::zeros);

  int K_star = K_init;

  arma::vec xi(K_max, fill::zeros);
  if(ind_slice){
    for(int k = 0; k < K_max; ++k)
      xi(k) = (1.0-rho) * std::pow(rho, (double)k);
  }

  // ========================== Initialisation ==========================
  if(init_A_pca){
    pca_init_X(X, d, U, A, psi);
  } else {
    A.zeros(p, d);
    for(int j = 0; j < p; ++j){
      double tau = R::rgamma(apsi, 1.0/bpsi);
      psi(j) = 1.0 / std::max(tau, 1e-12);
    }
  }
  for(int j = 0; j < p; ++j)
    if(!std::isfinite(psi(j)) || psi(j) <= 1e-12) psi(j) = 1e-6;

  for(int i = 0; i < n; ++i){
    double jit = std::fabs(R::rnorm(0.0, t_eps));
    t(i) = (Y(i)==1.0) ? (logC(i)-jit) : (logC(i)+jit);
  }

  if(init_z_kmeans) init_z_kmeans_helper(U, K_init, z);
  else              init_z_balanced_random_helper(n, K_init, z);

  arma::mat Psi0;
  arma::mat Psi0_base;
  if(init_A_pca){
    arma::mat U_cov = arma::cov(U);
    U_cov = 0.5 * (U_cov + U_cov.t()) +
            1e-6 * arma::eye<arma::mat>(d, d);
    double base_scale = (nu_w - d - 1.0) / std::max((double)K_init, 1.0);
    Psi0_base = base_scale * U_cov;
  } else {
    Psi0_base = arma::eye<arma::mat>(d, d);
  }

  if(auto_psi0){
    double base_tr = arma::trace(Psi0_base);
    if(base_tr <= 1e-12){
      Rcpp::warning("Psi0_base trace ~ 0; falling back to Psi0_scale=1.0");
      Psi0_scale = 1.0;
    } else {
      Psi0_scale = (double)d / base_tr;
    }
    Rcpp::Rcout << "[auto Psi0] base_trace=" << base_tr
                << "  -> Psi0_scale=" << Psi0_scale << "\n";
  }

  Psi0 = Psi0_scale * Psi0_base;
  Psi0 = 0.5 * (Psi0 + Psi0.t()) +
         1e-6 * arma::eye<arma::mat>(d, d);

  {
    arma::ivec nk(K_init, fill::zeros);
    for(int i = 0; i < n; ++i){
      int k = z(i);
      if(k >= 0 && k < K_init) nk(k)++;
    }
    int cum = 0;
    for(int k = 0; k < K_init; ++k){
      int mk = n - cum;
      double vv = R::rbeta(1.0+nk(k), std::max(alpha+(double)(mk-nk(k)), 1e-6));
      V(k) = clamp01(vv, 1e-12, 1.0-1e-12);
      cum += nk(k);
    }
    for(int k = K_init; k < K_max; ++k)
      V(k) = clamp01(R::rbeta(1.0, alpha), 1e-12, 1.0-1e-12);
  }

  {
    double prod = 1.0;
    for(int k = 0; k < K_max; ++k){
      w_vec(k) = std::max(V(k)*prod, 1e-300);
      prod *= (1.0-V(k));
    }
  }

  for(int i = 0; i < n; ++i)
    s(i) = ind_slice ? R::runif(0.0, xi(z(i)))
                     : R::runif(0.0, w_vec(z(i)));

  for(int k = 0; k < K_max; ++k){
    Sigma.slice(k) = Psi0 / std::max(nu_w - d - 1.0, 1.0);
    nu_mat.row(k)  = arma::zeros<arma::rowvec>(d);
    sig2(k) = 1.0/R::rgamma(a0, 1.0/b0);
    mu(k)   = R::rnorm(mu0, std::sqrt(sig2(k)/kappa_mu0));
  }
  for(int k = 0; k < K_init; ++k){
    arma::uvec idx = arma::find(z == k);
    if((int)idx.n_elem >= 2){
      arma::vec tk = t.elem(idx);
      mu(k)   = arma::mean(tk);
      double v = arma::as_scalar(arma::var(tk, 1));
      sig2(k) = (std::isfinite(v) && v > 1e-4) ? v : 1e-4;
    }
    if(!std::isfinite(sig2(k)) || sig2(k) <= 1e-12) sig2(k) = 1e-12;
  }

  delta.ones();

  // ========================== Output storage ==========================
  const int out_keep = (iters-burn)/thin;
  if(out_keep <= 0) stop("No kept draws: check (iters, burn, thin).");

  arma::mat  A_sum(p, d, fill::zeros);
  arma::vec  psi_sum(p, fill::zeros);
  arma::vec  delta_sum(p, fill::zeros);
  arma::mat  nu_sum(K_max, d, fill::zeros);
  arma::cube Sigma_sum(d, d, K_max, fill::zeros);
  arma::mat  U_sum(n, d, fill::zeros);
  arma::vec  beta_sum((q > 0 ? q : 1), fill::zeros);

  arma::mat      mu_draws(K_max, out_keep, fill::zeros);
  arma::mat      sig2_draws(K_max, out_keep, fill::zeros);
  arma::mat      psi_draws(p, out_keep, fill::zeros);
  arma::Mat<int> z_draws(n, out_keep, fill::zeros);
  arma::mat      t_draws(n, out_keep, fill::zeros);
  arma::vec      alpha_draws(out_keep, fill::zeros);
  arma::vec      a_alpha_draws(out_keep, fill::zeros);
  arma::vec      b_alpha_draws(out_keep, fill::zeros);
  arma::mat      w_draws(K_max, out_keep, fill::zeros);
  arma::vec      Kplus_draws(out_keep, fill::zeros);
  arma::vec      Kstar_draws(out_keep, fill::zeros);
  arma::mat      beta_draws((q > 0 ? q : 1), out_keep, fill::zeros);

  arma::cube nu_draws_c(K_max, d, out_keep, fill::zeros);
  arma::cube Sigma_draws_c(d, d, K_max*out_keep, fill::zeros);

  const int diag_keep_cap = std::min(out_keep, diag_max_keep);
  arma::cube U_draws_diag(n, d, diag_keep_cap, fill::zeros);
  arma::cube A_draws_diag(J, d, diag_keep_cap, fill::zeros);
  arma::mat  psi_draws_diag(J, diag_keep_cap, fill::zeros);

  int keep_idx = 0, keep_diag_idx = 0;

  // ========================== INIT print ==========================
  {
    arma::ivec nk0(K_init, fill::zeros);
    for(int i = 0; i < n; ++i){
      int k = z(i); if(k>=0 && k<K_init) nk0(k)++;
    }
    int Kp0 = 0; for(int k=0; k<K_init; ++k) if(nk0(k)>0) Kp0++;
    Rcpp::Rcout << "[INIT KGW v10 +ContSSVS] n=" << n << " p=" << p << " d=" << d
                << " q=" << q
                << " K_max=" << K_max << " K_init=" << K_init
                << " ind_slice=" << (ind_slice?1:0)
                << (ind_slice ? (" rho=" + std::to_string(rho)) : "")
                << " init_A_pca=" << (init_A_pca?1:0)
                << " init_z_kmeans=" << (init_z_kmeans?1:0)
                << " sample_alpha=" << (sample_alpha?1:0)
                << " alpha=" << alpha
                << " sample_ab=" << (sample_ab?1:0)
                << " a=" << a_alpha << " b=" << b_alpha << "\n";

    // NEW: SSVS mode banner
    Rcpp::Rcout << "[INIT] ssvs_mode=" << (use_continuous_ssvs ? "CONTINUOUS" : "POINT_MASS")
                << "  sigmaA2=" << sigmaA2;
    if(use_continuous_ssvs){
      Rcpp::Rcout << "  tau_spike2=" << tau_spike2
                  << "  (sigmaA/tau_spike=" << std::sqrt(sigmaA2/tau_spike2) << "x)";
    }
    Rcpp::Rcout << "  pi0=" << pi0 << "\n";

    if(q > 0){
      Rcpp::Rcout << "[INIT] beta initialized at prior mean b0_beta = [";
      for(int j = 0; j < q; ++j){
        Rcpp::Rcout << b0_beta(j);
        if(j < q-1) Rcpp::Rcout << ", ";
      }
      Rcpp::Rcout << "]\n";
    }
    Rcpp::Rcout << "[INIT] nu_w=" << nu_w
                << " Psi0_scale=" << Psi0_scale
                << (auto_psi0 ? " (AUTO)" : "")
                << " Psi0 trace=" << arma::trace(Psi0)
                << " (empirical Psi0 base=" << (init_A_pca?"YES (cov(U_pca))":"NO (identity)") << ")\n";
    Rcpp::Rcout << "[INIT] A: " << (init_A_pca?"PCA-initialized":"zero-initialized")
                << ",  psi median=" << arma::median(psi) << "\n";
    Rcpp::Rcout << "[INIT] Kplus=" << Kp0 << " w[1.." << K_init << "]:";
    for(int k=0; k<K_init; ++k) Rcpp::Rcout << " " << w_vec(k);
    Rcpp::Rcout << "\n";
  }

  const int PRINT_EVERY = 200;

  // ========================== Gibbs loop ==========================
  for(int it = 1; it <= iters; ++it){

    // (1) Slice variables
    if(ind_slice){
      for(int i = 0; i < n; ++i)
        s(i) = R::runif(0.0, xi(z(i)));
    } else {
      for(int i = 0; i < n; ++i)
        s(i) = R::runif(0.0, w_vec(z(i)));
    }

    // (2) Extend K*
    double s_star = arma::min(s);
    {
      if(ind_slice){
        double tail = 1.0;
        for(int k = 0; k < K_star; ++k) tail *= (1.0-V(k));
        while(K_star < K_max && xi(K_star) >= s_star){
          sig2(K_star) = 1.0/R::rgamma(a0, 1.0/b0);
          mu(K_star)   = R::rnorm(mu0, std::sqrt(sig2(K_star)/kappa_mu0));
          Sigma.slice(K_star) = Psi0 / std::max(nu_w - d - 1.0, 1.0);
          nu_mat.row(K_star)  = arma::zeros<arma::rowvec>(d);
          double vv = clamp01(R::rbeta(1.0, alpha), 1e-12, 1.0-1e-12);
          V(K_star) = vv;
          w_vec(K_star) = std::max(vv*tail, 1e-300);
          tail *= (1.0-vv);
          K_star++;
        }
      } else {
        double tail = 1.0;
        for(int k = 0; k < K_star; ++k) tail *= (1.0-V(k));
        while(tail >= s_star && K_star < K_max){
          double vv = clamp01(R::rbeta(1.0, alpha), 1e-12, 1.0-1e-12);
          V(K_star) = vv;
          w_vec(K_star) = std::max(vv*tail, 1e-300);
          tail *= (1.0-vv);
          sig2(K_star) = 1.0/R::rgamma(a0, 1.0/b0);
          mu(K_star)   = R::rnorm(mu0, std::sqrt(sig2(K_star)/kappa_mu0));
          Sigma.slice(K_star) = Psi0 / std::max(nu_w - d - 1.0, 1.0);
          nu_mat.row(K_star)  = arma::zeros<arma::rowvec>(d);
          K_star++;
        }
      }
    }

    // (3) Precompute Sigma_k^{-1}, log|Sigma_k|
    arma::cube Sigma_inv(d, d, K_star, fill::zeros);
    arma::vec  logdetS(K_star, fill::zeros);
    for(int k = 0; k < K_star; ++k){
      Sigma_inv.slice(k) = inv_sympd_ridge(Sigma.slice(k));
      logdetS(k)         = logdet_ridge(Sigma.slice(k));
    }

    arma::vec Wbeta(n, fill::zeros);
    if(q > 0) Wbeta = W * beta_cur;

    // (4) t_i
    for(int i = 0; i < n; ++i){
      int    k   = z(i);
      double m_i = Wbeta(i) + mu(k);
      double sd  = std::sqrt(std::max(sig2(k), 1e-12));
      double bnd = logC(i);
      t(i) = (Y(i)==1.0)
               ? rtruncnorm_right_icdf(m_i, sd, bnd)
               : rtruncnorm_left_icdf (m_i, sd, bnd);
    }

    // (5) z_i
    for(int i = 0; i < n; ++i){
      std::vector<int> active_i;
      active_i.reserve(K_star);
      if(ind_slice){
        for(int k = 0; k < K_star; ++k)
          if(xi(k) > s(i)) active_i.push_back(k);
      } else {
        for(int k = 0; k < K_star; ++k)
          if(w_vec(k) > s(i)) active_i.push_back(k);
      }

      if(active_i.empty()) continue;

      const int Kact = (int)active_i.size();
      arma::vec lp(Kact);
      for(int j = 0; j < Kact; ++j){
        int k = active_i[j];
        double m_tk = Wbeta(i) + mu(k);
        double lt = dnorm_log(t(i), m_tk, std::sqrt(std::max(sig2(k), 1e-12)));
        arma::rowvec diff = U.row(i) - nu_mat.row(k);
        double quad = arma::as_scalar(diff * Sigma_inv.slice(k) * diff.t());
        double lu   = -0.5*(d*std::log(2.0*arma::datum::pi) + logdetS(k) + quad);
        double lw;
        if(ind_slice){
          lw = std::log(std::max(w_vec(k), 1e-300))
               - std::log(xi(k));
        } else {
          lw = 0.0;
        }
        lp(j) = lw + lt + lu;
      }

      double mx = lp.max();
      arma::vec wt = arma::exp(lp-mx);
      double sw = arma::accu(wt);
      if(!std::isfinite(sw) || sw <= 0.0) continue;
      wt /= sw;
      double u = R::runif(0.0, 1.0), cum = 0.0;
      z(i) = active_i[Kact-1];
      for(int j = 0; j < Kact; ++j){
        cum += wt(j);
        if(u < cum){ z(i) = active_i[j]; break; }
      }
    }

    // (6) U_i
    for(int j = 0; j < p; ++j)
      if(!std::isfinite(psi(j)) || psi(j) <= 1e-12) psi(j) = 1e-6;

    arma::mat Psi_inv   = arma::diagmat(1.0/psi);
    arma::mat At_PsiInv = A.t() * Psi_inv;
    arma::mat B         = At_PsiInv * A;

    for(int i = 0; i < n; ++i){
      int       k    = z(i);
      arma::mat Prec = B + Sigma_inv.slice(k);
      arma::mat Cov  = inv_sympd_ridge(Prec);
      arma::vec rhs  = At_PsiInv * X.row(i).t()
                       + Sigma_inv.slice(k) * nu_mat.row(k).t();
      U.row(i) = rmvnorm(Cov*rhs, Cov).t();
    }

    // ===========================================================
    // (7) A, delta, psi | U, X     [CHANGED for continuous SSVS]
    // ===========================================================
    {
      arma::mat UtU = U.t() * U;
      arma::mat UtX = U.t() * X;

      // ---- (7a) delta_j update ----
      //
      // Marginal log-likelihood of x_j under prior A_{j.} ~ N(0, sigma_c^2 * I):
      //   log p(x_j | delta_j=c) = const - 0.5 log|G_c| + 0.5 (sigma_c^2 / psi_j^2) quad_c
      // where
      //   G_c    = I + (sigma_c^2 / psi_j) U'U
      //   quad_c = (U'x_j)' G_c^{-1} (U'x_j)
      //
      // POINT-MASS  (tau_spike2 <= 0):  sigma_0 = 0  ->  G_0 = I, log|G_0| = 0, quad_0 term vanishes
      // CONTINUOUS  (tau_spike2 > 0) :  sigma_0 = sqrt(tau_spike2), full formula
      //
      for(int j = 0; j < p; ++j){
        double psi_j = std::max(psi(j), 1e-12);
        arma::vec Utx = UtX.col(j);

        // Slab (delta=1) marginal terms
        double c1       = sigmaA2 / psi_j;
        arma::mat G1    = arma::eye<arma::mat>(d, d) + c1 * UtU;
        double logdetG1 = logdet_ridge(G1);
        arma::mat G1inv = inv_sympd_ridge(G1);
        double quad1    = arma::as_scalar(Utx.t() * G1inv * Utx);

        double logit;
        if(use_continuous_ssvs){
          // Spike (delta=0) marginal terms — same form, narrow prior
          double c0       = tau_spike2 / psi_j;
          arma::mat G0    = arma::eye<arma::mat>(d, d) + c0 * UtU;
          double logdetG0 = logdet_ridge(G0);
          arma::mat G0inv = inv_sympd_ridge(G0);
          double quad0    = arma::as_scalar(Utx.t() * G0inv * Utx);

          // log[ P(delta=1 | x_j) / P(delta=0 | x_j) ] =
          //   logit_pi0
          //   + 0.5 (log|G_0| - log|G_1|)
          //   + 0.5 / psi_j^2 * (sigmaA2 * quad_1 - tau_spike2 * quad_0)
          logit = logit_pi0
                + 0.5 * (logdetG0 - logdetG1)
                + 0.5 / (psi_j * psi_j) * (sigmaA2 * quad1 - tau_spike2 * quad0);
        } else {
          // Point-mass limit
          logit = logit_pi0 + 0.5 * (c1 / psi_j * quad1 - logdetG1);
        }
        if(!std::isfinite(logit)) logit = (logit > 0) ? 30.0 : -30.0;
        double pr1 = 1.0 / (1.0 + std::exp(-logit));
        delta(j) = (R::runif(0.0, 1.0) < pr1) ? 1.0 : 0.0;
      }

      // ---- (7b) A_j update ----
      //
      // POINT-MASS:  delta_j = 0 -> A_{j.} = 0 ; delta_j = 1 -> draw from slab posterior
      // CONTINUOUS:  always draw, with prior variance = sigmaA2 (delta=1) or tau_spike2 (delta=0)
      //
      for(int j = 0; j < p; ++j){
        double psi_j     = std::max(psi(j), 1e-12);
        double prior_var;

        if(use_continuous_ssvs){
          prior_var = (delta(j) >= 0.5) ? sigmaA2 : tau_spike2;
        } else {
          if(delta(j) < 0.5){ A.row(j).zeros(); continue; }
          prior_var = sigmaA2;
        }

        arma::mat Prec = (1.0 / psi_j) * UtU
                       + (1.0 / prior_var) * arma::eye<arma::mat>(d, d);
        arma::mat Cov  = inv_sympd_ridge(Prec);
        arma::vec mean = Cov * ((1.0 / psi_j) * UtX.col(j));
        A.row(j)       = rmvnorm(mean, Cov).t();
      }

      // ---- (7c) psi_j update (unchanged) ----
      for(int j = 0; j < p; ++j){
        arma::vec r  = X.col(j) - U * A.row(j).t();
        double    ss = arma::dot(r, r);
        double    tau = R::rgamma(apsi + 0.5 * n, 1.0 / (bpsi + 0.5 * ss));
        psi(j) = 1.0 / std::max(tau, 1e-12);
        if(!std::isfinite(psi(j)) || psi(j) <= 0) psi(j) = 1e-6;
      }
    }

    // (8) mu_k, sig2_k
    for(int k = 0; k < K_star; ++k){
      arma::uvec idx = arma::find(z == k);
      int nk = (int)idx.n_elem;
      if(nk == 0){
        sig2(k) = 1.0/R::rgamma(a0, 1.0/b0);
        mu(k)   = R::rnorm(mu0, std::sqrt(sig2(k)/kappa_mu0));
      } else {
        arma::vec r_k = t.elem(idx) - Wbeta.elem(idx);
        double rbar   = arma::mean(r_k);
        double sse    = arma::accu(arma::square(r_k-rbar));
        double add    = (kappa_mu0*nk)/(kappa_mu0+nk)*std::pow(rbar-mu0,2.0);
        sig2(k) = 1.0/std::max(R::rgamma(a0+0.5*nk, 1.0/(b0+0.5*sse+0.5*add)), 1e-12);
        double kn = kappa_mu0+nk;
        mu(k) = R::rnorm((kappa_mu0*mu0+nk*rbar)/kn, std::sqrt(sig2(k)/kn));
      }
      if(!std::isfinite(sig2(k)) || sig2(k)<=1e-12) sig2(k) = 1e-12;
    }

    // (9) beta
    if(q > 0){
      arma::vec r_vec(n), wts(n);
      for(int i = 0; i < n; ++i){
        r_vec(i) = t(i) - mu(z(i));
        wts(i)   = 1.0 / std::max(sig2(z(i)), 1e-12);
      }
      arma::mat Wtil = W.each_col() % arma::sqrt(wts);
      arma::mat WDinvW = Wtil.t() * Wtil;
      arma::vec WDinvr = W.t() * (wts % r_vec);
      arma::mat B_n    = inv_sympd_ridge(B0_inv + WDinvW);
      arma::vec b_n    = B_n * (B0inv_b0 + WDinvr);
      B_n = 0.5 * (B_n + B_n.t());
      beta_cur = rmvnorm(b_n, B_n);
    }

    // (10) V_k, w_k, alpha
    {
      arma::ivec nk(K_star, fill::zeros);
      for(int i = 0; i < n; ++i){
        int k = z(i);
        if(k >= 0 && k < K_star) nk(k)++;
      }

      int cum = 0;
      for(int k = 0; k < K_star; ++k){
        int mk = n - cum;
        double b_par = std::max(alpha + (double)(mk - nk(k)), 1e-6);
        double vv = R::rbeta(1.0 + nk(k), b_par);
        V(k) = clamp01(vv, 1e-12, 1.0-1e-12);
        cum += nk(k);
      }

      double prod = 1.0;
      for(int k = 0; k < K_star; ++k){
        w_vec(k) = std::max(V(k)*prod, 1e-300);
        prod *= (1.0-V(k));
      }
      for(int k = K_star; k < K_max; ++k){
        double vv = clamp01(R::rbeta(1.0, alpha), 1e-12, 1.0-1e-12);
        V(k) = vv;
        w_vec(k) = std::max(vv*prod, 1e-300);
        prod *= (1.0-vv);
      }
      {
        int K_new = 1;
        for(int i = 0; i < n; ++i){
          int k = z(i);
          if(k + 1 > K_new) K_new = k + 1;
        }
        K_star = K_new;
      }

      if(sample_alpha){
        arma::ivec nk2(K_star, fill::zeros);
        for(int i = 0; i < n; ++i){ int k=z(i); if(k>=0&&k<K_star) nk2(k)++; }
        int K_plus = 0;
        for(int k = 0; k < K_star; ++k) if(nk2(k)>0) K_plus++;
        update_alpha_ew(alpha, n, K_plus, a_alpha, b_alpha);

        if(sample_ab){
          update_a_hyperprior(a_alpha, alpha, b_alpha,
                              alpha_a_prior, beta_a_prior,
                              mh_a_sd, a_accept_ct, a_total_ct);
          update_b_hyperprior(b_alpha, alpha, a_alpha,
                              alpha_b_prior, beta_b_prior);
        }
      }
    }

    // (11) nu_k, Sigma_k
    for(int k = 0; k < K_star; ++k){
      arma::uvec idx = arma::find(z == k);
      int nk = (int)idx.n_elem;
      if(nk == 0){
        Sigma.slice(k) = rinvwishart(Psi0, (int)std::round(nu_w));
        nu_mat.row(k)  = rmvnorm(nu0, Sigma.slice(k)/kappa0).t();
        continue;
      }
      arma::mat Uk      = U.rows(idx);
      arma::rowvec ubar = arma::mean(Uk, 0);
      arma::mat Sk(d, d, fill::zeros);
      for(int ii = 0; ii < nk; ++ii){
        arma::vec diff = Uk.row(ii).t()-ubar.t();
        Sk += diff*diff.t();
      }
      arma::vec dm    = ubar.t()-nu0;
      arma::mat Psi_n = Psi0+Sk+(kappa0*nk/(kappa0+nk))*(dm*dm.t());
      Sigma.slice(k) = rinvwishart(Psi_n, (int)std::round(nu_w+nk));
      nu_mat.row(k)  = rmvnorm((kappa0*nu0+nk*ubar.t())/(kappa0+nk),
                                Sigma.slice(k)/(kappa0+nk)).t();
    }

    // ---- progress print ----
    if(it==1 || it%PRINT_EVERY==0 || it==burn || it==iters){
      arma::ivec nk(K_star, fill::zeros);
      for(int i=0; i<n; ++i){ int k=z(i); if(k>=0&&k<K_star) nk(k)++; }
      int Kplus=0; for(int k=0; k<K_star; ++k) if(nk(k)>0) Kplus++;
      Rcpp::Rcout << "[it " << it << "] K*=" << K_star
                  << " Kplus=" << Kplus
                  << " s*=" << s_star
                  << " alpha=" << alpha;
      if(sample_ab){
        double acc = (a_total_ct > 0) ? (double)a_accept_ct/(double)a_total_ct : 0.0;
        Rcpp::Rcout << " a=" << a_alpha << " b=" << b_alpha
                    << " acc_a=" << std::fixed << std::setprecision(2) << acc;
        Rcpp::Rcout.unsetf(std::ios_base::floatfield);
      }
      if(q > 0){
        Rcpp::Rcout << " beta=[";
        Rcpp::Rcout << std::fixed << std::setprecision(3);
        for(int j = 0; j < q; ++j){
          Rcpp::Rcout << beta_cur(j);
          if(j < q-1) Rcpp::Rcout << ",";
        }
        Rcpp::Rcout << "]";
        Rcpp::Rcout.unsetf(std::ios_base::floatfield);
      }
      int n_active = 0;
      for(int j = 0; j < p; ++j) if(delta(j) > 0.5) n_active++;
      Rcpp::Rcout << " #delta=1:" << n_active;
      Rcpp::Rcout << " n_k:";
      for(int k=0; k<K_star; ++k) Rcpp::Rcout << " " << nk(k);
      Rcpp::Rcout << "\n";
    }

    // ---- store kept draws ----
    if(it > burn && ((it-burn)%thin == 0)){
      if(keep_idx >= out_keep) break;

      A_sum     += A;
      psi_sum   += psi;
      delta_sum += delta;
      for(int k=0; k<K_star; ++k){
        nu_sum.row(k)      += nu_mat.row(k);
        Sigma_sum.slice(k) += Sigma.slice(k);
      }
      U_sum += U;
      if(q > 0) beta_sum += beta_cur;

      mu_draws.col(keep_idx).subvec(0,K_star-1)   = mu.subvec(0,K_star-1);
      sig2_draws.col(keep_idx).subvec(0,K_star-1) = sig2.subvec(0,K_star-1);
      psi_draws.col(keep_idx)   = psi;
      z_draws.col(keep_idx)     = z;
      t_draws.col(keep_idx)     = t;
      alpha_draws(keep_idx)     = alpha;
      a_alpha_draws(keep_idx)   = a_alpha;
      b_alpha_draws(keep_idx)   = b_alpha;
      w_draws.col(keep_idx).subvec(0,K_star-1) = w_vec.subvec(0,K_star-1);
      if(q > 0) beta_draws.col(keep_idx) = beta_cur;

      nu_draws_c.slice(keep_idx) = nu_mat;
      for(int k = 0; k < K_max; ++k)
        Sigma_draws_c.slice(k*out_keep + keep_idx) = Sigma.slice(k);

      arma::ivec nk(K_star, fill::zeros);
      for(int i=0; i<n; ++i){ int k=z(i); if(k>=0&&k<K_star) nk(k)++; }
      int Kplus=0; for(int k=0; k<K_star; ++k) if(nk(k)>0) Kplus++;

      Kplus_draws(keep_idx) = (double)Kplus;
      Kstar_draws(keep_idx) = (double)K_star;

      if(keep_diag_idx < diag_keep_cap){
        U_draws_diag.slice(keep_diag_idx) = U;
        for(int jj=0; jj<J; ++jj){
          int j0 = diag0[jj];
          A_draws_diag.slice(keep_diag_idx).row(jj) = A.row(j0);
          psi_draws_diag(jj,keep_diag_idx) = psi(j0);
        }
        keep_diag_idx++;
      }
      keep_idx++;
    }
  }

  if(keep_idx <= 0) stop("No kept draws — check (iters, burn, thin).");

  arma::Mat<int> z_out = z_draws.cols(0, keep_idx-1);
  z_out += 1;

  IntegerVector dfeat(J);
  for(int jj=0; jj<J; ++jj) dfeat[jj] = diag0[jj]+1;

  List out;
  out["A_mean"]      = A_sum   / (double)keep_idx;
  out["psi_mean"]    = psi_sum / (double)keep_idx;
  out["pip"]         = delta_sum / (double)keep_idx;
  out["U_mean"]      = U_sum   / (double)keep_idx;
  out["nu_mean"]     = nu_sum  / (double)keep_idx;
  out["Sigma_mean"]  = Sigma_sum / (double)keep_idx;
  out["kept"]        = keep_idx;

  out["mu_draws"]    = mu_draws.cols(0, keep_idx-1);
  out["sig2_draws"]  = sig2_draws.cols(0, keep_idx-1);
  out["psi_draws"]   = psi_draws.cols(0, keep_idx-1);
  out["z_draws"]     = z_out;
  out["t_draws"]     = t_draws.cols(0, keep_idx-1);
  out["alpha_draws"] = alpha_draws.subvec(0, keep_idx-1);
  out["a_alpha_draws"] = a_alpha_draws.subvec(0, keep_idx-1);
  out["b_alpha_draws"] = b_alpha_draws.subvec(0, keep_idx-1);
  out["a_alpha_accept_rate"] = (a_total_ct > 0) ? (double)a_accept_ct/(double)a_total_ct : 0.0;
  out["w_draws"]     = w_draws.cols(0, keep_idx-1);
  out["Kplus_draws"] = Kplus_draws.subvec(0, keep_idx-1);
  out["Kstar_draws"] = Kstar_draws.subvec(0, keep_idx-1);

  {
    arma::cube nu_out(K_max, d, keep_idx, fill::zeros);
    for(int s = 0; s < keep_idx; ++s)
      nu_out.slice(s) = nu_draws_c.slice(s);
    out["nu_draws"] = nu_out;
  }
  {
    arma::cube Sigma_out(d, d, K_max*keep_idx, fill::zeros);
    for(int k = 0; k < K_max; ++k)
      for(int s = 0; s < keep_idx; ++s)
        Sigma_out.slice(k*keep_idx + s) = Sigma_draws_c.slice(k*out_keep + s);
    out["Sigma_draws"] = Sigma_out;
  }

  out["diag_feat_idx"]   = dfeat;
  out["diag_kept"]       = keep_diag_idx;
  out["U_draws_diag"]    = U_draws_diag.slices(0, keep_diag_idx-1);
  out["A_draws_diag"]    = A_draws_diag.slices(0, keep_diag_idx-1);
  out["psi_draws_diag"]  = psi_draws_diag.cols(0, keep_diag_idx-1);

  out["ind_slice"]       = ind_slice;
  out["rho"]             = rho;
  out["K_max"]           = K_max;
  out["K_init"]          = K_init;
  out["init_A_pca"]      = init_A_pca;
  out["init_z_kmeans"]   = init_z_kmeans;
  out["init_z"]          = std::string(init_z_kmeans ? "kmeans_pca" : "balanced_random");

  out["q"]               = q;
  out["Psi0"]            = Psi0;
  out["empirical_Psi0"]  = init_A_pca;
  out["nu_w"]            = nu_w;
  out["Psi0_scale"]      = Psi0_scale;

  // NEW: continuous SSVS provenance
  out["ssvs_mode"]       = std::string(use_continuous_ssvs ? "continuous" : "point_mass");
  out["sigmaA2"]         = sigmaA2;
  out["tau_spike2"]      = tau_spike2;
  if(use_continuous_ssvs){
    out["spike_slab_sd_ratio"] = std::sqrt(sigmaA2 / std::max(tau_spike2, 1e-12));
  }

  if(q > 0){
    out["beta_mean"]  = beta_sum / (double)keep_idx;
    out["beta_draws"] = beta_draws.cols(0, keep_idx-1);
  }

  return out;
}
