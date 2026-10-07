# Preprocess x, u, and fixed a
prepare_data <- function(data) {
  if (is.null(data$x) || is.null(data$u) || is.null(data$a)) {
    stop("data must contain x, u, and a.");
  }
  
  x <- as.matrix(data$x);
  u <- as.matrix(data$u);
  a <- as.matrix(data$a);
  storage.mode(x) <- "double";
  storage.mode(u) <- "double";
  storage.mode(a) <- "double";
  
  # Cells, genes, regulators
  N <- nrow(x);
  J <- ncol(x);
  K <- if (is.null(data$K)) ncol(a) else as.integer(data$K);
  
  data$x <- x;
  data$u <- u;
  data$a <- a;
  data$N <- N;
  data$J <- J;
  data$K <- K;
  
  # C_ik = (1-u_ik) a_ik
  data$C <- (1 - u[, seq_len(K), drop=FALSE]) * a;
  
  structure(data, class="exprsim_vi_data")
}

# Initialize variational parameters
prepare_params <- function(params0, data, hparams) {
  if (is.null(params0)) {
    params0 <- list();
  }
  
  N <- data$N;
  J <- data$J;
  K <- data$K;
  
  # Start Beta variational factors at their priors
  vparams <- list(
    psi = list(a = rep(hparams$alpha_psi, K), b = rep(hparams$beta_psi, K)),
    phi = list(a = rep(hparams$alpha_phi, J), b = rep(hparams$beta_phi, J)),
    eta = list(a = rep(hparams$alpha_eta, J), b = rep(hparams$beta_eta, J)),
    beta_delta = list(
      mu_beta = if (is.null(params0$mu_beta)) {
        matrix(0, K, J)
      } else {
        as.matrix(params0$mu_beta)
      },
      sigma2_beta = if (is.null(params0$sigma2_beta)) {
        matrix(1 / hparams$tau, K, J)
      } else {
        as.matrix(params0$sigma2_beta)
      }
    ),
    z = list(a = matrix(1, N, J), b = matrix(hparams$gamma, N, J))
  );
  
  # Initial E[delta_kj] = E[psi_k] E[phi_j]
  if (is.null(params0$expect_delta)) {
    expect_psi <- expect_x_beta(vparams$psi);
    expect_phi <- expect_x_beta(vparams$phi);
    expect_delta <- outer(expect_psi, expect_phi, "*");
  } else {
    expect_delta <- as.matrix(params0$expect_delta);
  }
  
  # Initial E[v_ij] = E[eta_j], with v_ij = 1 when x_ij > 0
  if (is.null(params0$expect_v)) {
    expect_eta <- expect_x_beta(vparams$eta);
    expect_v <- matrix(rep(expect_eta, each=N), N, J);
    expect_v[data$x > 0] <- 1;
  } else {
    expect_v <- as.matrix(params0$expect_v);
    expect_v[data$x > 0] <- 1;
  }
  
  # Initial auxiliary tangent parameter for the z lower bound
  if (is.null(params0$xi_z)) {
    xi_z <- matrix(hparams$rho, N, J) + (1 - data$u) / hparams$gamma;
  } else {
    xi_z <- as.matrix(params0$xi_z);
  }
  
  # Initial auxiliary parameter for the delta Jensen lower bound
  if (is.null(params0$xi_delta)) {
    expect_log_psi <- expect_log_x_beta(vparams$psi);
    expect_log_one_minus_psi <- expect_log_1mx_beta(vparams$psi);
    expect_log_one_minus_phi <- expect_log_1mx_beta(vparams$phi);
    xi_delta <- logistic(outer(expect_log_one_minus_psi - expect_log_psi, expect_log_one_minus_phi, "-"));
  } else {
    xi_delta <- as.matrix(params0$xi_delta);
  }
  
  mu_beta <- vparams$beta_delta$mu_beta;
  sigma2_beta <- vparams$beta_delta$sigma2_beta;
  
  params <- list(
    expect_delta = expect_delta,
    expect_delta_beta = expect_delta * mu_beta,
    expect_beta = expect_delta * mu_beta,
    expect_beta2 = (1 - expect_delta) / hparams$tau + expect_delta * (mu_beta^2 + sigma2_beta),
    expect_v = expect_v,
    xi_z = xi_z,
    xi_delta = xi_delta
  );
  
  list(params = params, vparams = vparams)
}

# Decide which variational blocks are fixed during VI
prepare_fixed <- function(fixed) {
  if (is.null(fixed)) {
    fixed <- list();
  }
  
  for (block_name in c("beta_delta", "z_v", "eta", "psi", "phi")) {
    if (is.null(fixed[[block_name]])) {
      fixed[[block_name]] <- FALSE;
    }
  }
  
  # First implementation: a is fixed
  if (is.null(fixed$a)) {
    fixed$a <- TRUE;
  }
  
  fixed
}