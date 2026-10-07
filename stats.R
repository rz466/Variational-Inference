### Functions from github that are reused directly

# E[x] under a beta distribution
expect_x_beta <- function(param) {
  with(param, a / (a + b))
}

# E[logit(x)] under a beta distribution
expect_logit_x_beta <- function(param) {
  with(param, digamma(a) - digamma(b))
}

# E[log(x)] under a beta distribution
expect_log_x_beta <- function(param) {
  with(param, digamma(a) - digamma(a + b))
}

# E[log(1-x)] under a beta distribution
expect_log_1mx_beta <- function(param) {
  with(param, digamma(b) - digamma(a + b))
}

# -E[log q] for a Bernoulli distribution, with direct parameterization
entropy_bern_d <- function(prob) {
  # -prob*log(prob) - (1-prob)*log(1-prob)
  log_q_zero <- ifelse(prob == 1, 0, -(1 - prob) * log(1 - prob));
  log_q_one <- ifelse(prob == 0, 0, -prob * log(prob));
  sum(log_q_zero + log_q_one)
}

### Newly added

# Elementwise version of -E[log q] for a Bernoulli distribution
bern_h <- function(prob) {
  out <- numeric(length(prob));
  idx_positive <- prob > 0;
  idx_less_than_one <- prob < 1;
  out[idx_positive] <- out[idx_positive] -
    prob[idx_positive] * log(prob[idx_positive]);
  out[idx_less_than_one] <- out[idx_less_than_one] -
    (1 - prob[idx_less_than_one]) * log(1 - prob[idx_less_than_one]);
  dim(out) <- dim(prob);
  out
}

# multiplier * log(log_argument), with 0*log(log_argument)=0
x_log_y <- function(multiplier, log_argument) {
  stopifnot(length(multiplier) == length(log_argument));
  out <- numeric(length(multiplier));
  idx_nonzero <- multiplier != 0;
  out[idx_nonzero] <- multiplier[idx_nonzero] * log(log_argument[idx_nonzero]);
  dim(out) <- dim(multiplier);
  out
}

# log E[exp(delta * exponent)] = log[(1-prob) + prob*exp(exponent)]
log_expect_exp_bernoulli <- function(prob, exponent) {
  if (prob <= 0) {
    return(rep(0, length(exponent)));
  }
  if (prob >= 1) {
    return(exponent);
  }
  log_delta_zero <- log1p(-prob);
  log_delta_one <- log(prob) + exponent;
  max_log_component <- pmax(log_delta_zero, log_delta_one);
  max_log_component + log(exp(log_delta_zero - max_log_component) + exp(log_delta_one - max_log_component))
}

# E(m), E(log(m)), E(m^-1), E(m*log(m))
compute_m_moments <- function(data, expect_delta, mu_beta, sigma2_beta) {
  N <- data$N;
  J <- data$J;
  K <- data$K;
  
  # C_ik = (1-u_ik) a_ik
  C_ik <- data$C;
  C_ik_squared <- C_ik * C_ik;
  
  expect_m <- matrix(0, N, J);
  expect_log_m <- matrix(0, N, J);
  expect_m_inverse <- matrix(0, N, J);
  expect_m_log_m <- matrix(0, N, J);
  
  for (j in seq_len(J)) {
    expect_delta_j <- expect_delta[, j];
    mu_beta_j <- mu_beta[, j];
    sigma2_beta_j <- sigma2_beta[, j];
    
    # E[m_ij]
    positive_exponent <- sweep(C_ik, 2, mu_beta_j, "*") + 0.5 * sweep(C_ik_squared, 2, sigma2_beta_j, "*");
    log_expect_exp_delta_beta <- matrix(0, N, K);
    
    for (k in seq_len(K)) {
      log_expect_exp_delta_beta[, k] <- log_expect_exp_bernoulli(expect_delta_j[k], positive_exponent[, k]);
    }
    
    log_expect_m <- rowSums(log_expect_exp_delta_beta);
    expect_m[, j] <- exp(log_expect_m);
    
    # E[log m_ij]
    expect_log_m[, j] <- as.vector(C_ik %*% (expect_delta_j * mu_beta_j));
    
    # E[m_ij log m_ij]
    expect_m_log_m_sum <- numeric(N);
    
    for (k in seq_len(K)) {
      if (expect_delta_j[k] > 0) {
        exp_beta_over_mixture <- exp(positive_exponent[, k] - log_expect_exp_delta_beta[, k]);
        expect_m_log_m_sum <- expect_m_log_m_sum + expect_delta_j[k] * C_ik[, k] *
          (mu_beta_j[k] + C_ik[, k] * sigma2_beta_j[k]) * exp_beta_over_mixture;
      }
    }
    
    expect_m_log_m[, j] <- expect_m[, j] * expect_m_log_m_sum;
    
    # E[m_ij^{-1}]
    negative_exponent <- sweep(C_ik, 2, -mu_beta_j, "*") + 0.5 * sweep(C_ik_squared, 2, sigma2_beta_j, "*");
    log_expect_exp_negative_delta_beta <- matrix(0, N, K);
    
    for (k in seq_len(K)) {
      log_expect_exp_negative_delta_beta[, k] <- log_expect_exp_bernoulli(expect_delta_j[k], negative_exponent[, k]);
    }
    
    expect_m_inverse[, j] <- exp(rowSums(log_expect_exp_negative_delta_beta));
  }
  
  if (
    any(!is.finite(expect_m)) ||
    any(!is.finite(expect_log_m)) ||
    any(!is.finite(expect_m_inverse)) ||
    any(!is.finite(expect_m_log_m))
  ) {
    stop(
      "Non-finite m-moment."
    );
  }
  
  list(
    expect_m = expect_m,
    expect_log_m = expect_log_m,
    expect_m_inverse = expect_m_inverse,
    expect_m_log_m = expect_m_log_m
  )
}

# E[log p(theta)] - E[log q(theta)] for a Beta prior and Beta variational distribution
beta_prior_minus_q <- function(variational_beta_param, prior_alpha, prior_beta) {
  expect_log_theta <- expect_log_x_beta(variational_beta_param);
  expect_log_one_minus_theta <- expect_log_1mx_beta(variational_beta_param);
  
  with(variational_beta_param, sum(-lbeta(prior_alpha, prior_beta) + lbeta(a, b) + (prior_alpha - a) * 
                                     expect_log_theta + (prior_beta - b) * expect_log_one_minus_theta))
}

# log p(a)
log_prior_a <- function(a, alpha) {
  sum(-lgamma(alpha) + alpha * log(alpha) + (alpha - 1) * log(a) - alpha * a)
}