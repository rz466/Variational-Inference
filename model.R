# Find a local maximum for the fixed-a ExprSim variational model
exprsim_vi_local <- function(data=NULL, hparams=NULL, params0=NULL, fixed=NULL,
                             model=NULL, verbose=1, max.iter=100, tol=0.001) {
  if (is.null(data)) {
    if (is.null(model)) {
      stop("Input data must be provided unless a previously fitted model is given.");
    }
  }
  
  if (is.null(model)) {
    model <- make_model(data, params0, hparams);
  } else {
    stopifnot(is(model, "exprsim_vi_fit"));
  }
  
  model$fixed <- prepare_fixed(fixed);
  
  if (!model$fixed$a) {
    stop("q(a) is not implemented in the fixed-a ExprSim VI model.");
  }
  
  niters <- 0;
  previous_elbo <- elbo(model);
  
  # ELBO change
  delta_elbo <- Inf;
  elbos <- rep(NA, max.iter);
  
  # Record the initial ELBO and its components
  initial_components <- elbo_components(model, include_fixed_a_prior=TRUE);
  
  model$elbo_trace <- rbind(
    model$elbo_trace,
    data.frame(
      iteration = 0,
      step = "initial",
      elbo = initial_components$total,
      delta = NA,
      z = initial_components$z,
      v_eta = initial_components$v_eta,
      beta_delta = initial_components$beta_delta,
      psi = initial_components$psi,
      phi = initial_components$phi,
      eta = initial_components$eta,
      a_fixed = initial_components$a_fixed,
      stringsAsFactors = FALSE
    )
  );
  
  # Convergence and stopping rule
  while ((!is.finite(delta_elbo) || abs(delta_elbo) > tol) && niters < max.iter) {
    niters <- niters + 1;
    step_elbo <- previous_elbo;
    
    # q(beta,delta): Taylor approximation might break monotonicity!
    if (!model$fixed$beta_delta) {
      model <- update_q_beta_delta(model);
      trace_result <- append_elbo_trace(
        model,
        niters,
        "q(beta,delta)",
        step_elbo,
        allow_decrease = TRUE,
        verbose = verbose
      );
      model <- trace_result$model;
      step_elbo <- trace_result$elbo;
    }
    
    # q(z,v)
    if (!model$fixed$z_v) {
      model <- update_q_z_v(model);
      trace_result <- append_elbo_trace(
        model,
        niters,
        "q(z,v)",
        step_elbo,
        verbose = verbose
      );
      model <- trace_result$model;
      step_elbo <- trace_result$elbo;
    }
    
    # xi_ij
    model <- update_xi_z(model);
    trace_result <- append_elbo_trace(
      model,
      niters,
      "xi_z",
      step_elbo,
      verbose = verbose
    );
    model <- trace_result$model;
    step_elbo <- trace_result$elbo;
    
    # q(eta)
    if (!model$fixed$eta) {
      model <- update_q_eta(model);
      trace_result <- append_elbo_trace(
        model,
        niters,
        "q(eta)",
        step_elbo,
        verbose = verbose
      );
      model <- trace_result$model;
      step_elbo <- trace_result$elbo;
    }
    
    # q(psi)
    if (!model$fixed$psi) {
      model <- update_q_psi(model);
      trace_result <- append_elbo_trace(
        model,
        niters,
        "q(psi)",
        step_elbo,
        verbose = verbose
      );
      model <- trace_result$model;
      step_elbo <- trace_result$elbo;
    }
    
    # xi_kj after q(psi)
    model <- update_xi_delta(model);
    trace_result <- append_elbo_trace(
      model,
      niters,
      "xi_delta_after_psi",
      step_elbo,
      verbose = verbose
    );
    model <- trace_result$model;
    step_elbo <- trace_result$elbo;
    
    # q(phi)
    if (!model$fixed$phi) {
      model <- update_q_phi(model);
      trace_result <- append_elbo_trace(
        model,
        niters,
        "q(phi)",
        step_elbo,
        verbose = verbose
      );
      model <- trace_result$model;
      step_elbo <- trace_result$elbo;
    }
    
    # xi_kj after q(phi)
    model <- update_xi_delta(model);
    trace_result <- append_elbo_trace(
      model,
      niters,
      "xi_delta_after_phi",
      step_elbo,
      verbose = verbose
    );
    model <- trace_result$model;
    step_elbo <- trace_result$elbo;
    
    # End of one complete VI iteration
    current_elbo <- step_elbo;
    delta_elbo <- current_elbo - previous_elbo;
    elbos[niters] <- current_elbo;
    
    if (verbose >= 3) {
      message("update iteration: ", niters, ",  elbo: ", current_elbo);
      
      if (verbose >= 4) {
        message("delta elbo: ", delta_elbo);
      }
    }
    
    previous_elbo <- current_elbo;
  }
  
  if (verbose >= 2) {
    message("elapsed iterations: ", niters);
  }
  
  if (niters > 0) {
    model$elbo <- c(model$elbo, elbos[seq_len(niters)]);
  }
  
  model
}

# Update E[m_ij], E[log m_ij], E[m_ij^-1], and E[m_ij log m_ij]
update_m_moments <- function(model) {
  m_moments <- compute_m_moments(
    data = model$data,
    expect_delta = model$params$expect_delta,
    mu_beta = model$vparams$beta_delta$mu_beta,
    sigma2_beta = model$vparams$beta_delta$sigma2_beta
  );
  
  model$params$expect_m <- m_moments$expect_m;
  model$params$expect_log_m <- m_moments$expect_log_m;
  model$params$expect_m_inverse <- m_moments$expect_m_inverse;
  model$params$expect_m_log_m <- m_moments$expect_m_log_m;
  
  model
}

# Update the Jensen auxiliary parameter xi_kj
update_xi_delta <- function(model) {
  model$params$xi_delta <- logistic(
    outer(model$params$expect_log_one_minus_psi - model$params$expect_log_psi, 
          model$params$expect_log_one_minus_phi, "-")
  );
  
  model
}

# Update the tangent auxiliary parameter xi_ij
update_xi_z <- function(model) {
  model$params$xi_z <- model$params$rho + (1 - model$data$u) * (model$vparams$z$a / model$vparams$z$b);
  
  model
}

# Update q(beta_kj, delta_kj)
update_q_beta_delta <- function(model) {
  J <- model$data$J;
  K <- model$data$K;
  C_ik <- model$data$C;
  C_ik_squared <- C_ik * C_ik;
  gamma <- model$params$gamma;
  tau <- model$params$tau;
  
  for (j in seq_len(J)) {
    
    # E[log z_ij | v_ij=1]
    expect_log_z_given_v_one <- digamma(model$vparams$z$a[, j]) - log(model$vparams$z$b[, j]);
    
    # Loop through perturbations
    for (k in seq_len(K)) {
      
      # Current E[n_ij] = E[log m_ij]
      expect_n <- model$params$expect_log_m[, j];
      gamma_exp_expect_n <- gamma * exp(expect_n);
      
      if (any(!is.finite(gamma_exp_expect_n)) || any(gamma_exp_expect_n <= 0)) {
        stop("Non-finite gamma*exp(E[n]) in q(beta,delta) update.");
      }
      
      # First derivative
      f_first <- model$params$expect_v[, j] * gamma_exp_expect_n * (log(gamma) + expect_log_z_given_v_one - digamma(gamma_exp_expect_n));
      
      # Second derivative
      f_second <- model$params$expect_v[, j] * gamma_exp_expect_n * (log(gamma) + expect_log_z_given_v_one - digamma(gamma_exp_expect_n) -
                                                                       gamma_exp_expect_n * trigamma(gamma_exp_expect_n));
      
      # Save old E[delta_kj beta_kj]
      old_expect_delta_beta <- model$params$expect_delta_beta[k, j];
      C_k <- C_ik[, k];
      C_k_squared <- C_ik_squared[, k];
      
      # delta_{beta_kj}^2
      precision_beta <- tau - sum(C_k_squared * f_second);
      
      if (!is.finite(precision_beta) || precision_beta <= 0) {
        stop("Taylor q(beta,delta) update is not normalizable at k=", k, ", j=", j, ": precision = ", precision_beta, ".");
      }
      
      sigma2_beta <- 1 / precision_beta;
      
      # mu_{beta_kj}
      mu_beta <- sigma2_beta * sum(C_k * f_first - C_k_squared * old_expect_delta_beta * f_second);
      
      # Jensen auxiliary parameter
      xi_delta <- model$params$xi_delta[k, j];
      H_xi_delta <- bern_h(xi_delta);
      
      # delta_kj = 0: log c_0
      log_c_zero <-
        xi_delta * model$params$expect_log_one_minus_psi[k] +
        (1 - xi_delta) * model$params$expect_log_psi[k] +
        (1 - xi_delta) * model$params$expect_log_one_minus_phi[j] +
        H_xi_delta + 0.5 * (log(2 * pi) - log(tau));
      
      # delta_kj = 1: log c_1
      log_c_one <-
        model$params$expect_log_psi[k] +
        model$params$expect_log_phi[j] +
        mu_beta^2 / (2 * sigma2_beta) +
        0.5 * log(2 * pi * sigma2_beta);
      
      # Normalize the two weights
      log_denominator <- log_sum_exp(c(log_c_zero, log_c_one));
      
      # E[delta_kj]
      expect_delta <- exp(log_c_one - log_denominator);
      
      # Save the new q(beta,delta) parameters
      model$vparams$beta_delta$mu_beta[k, j] <- mu_beta;
      model$vparams$beta_delta$sigma2_beta[k, j] <- sigma2_beta;
      model$params$expect_delta[k, j] <- expect_delta;
      
      # E[delta_kj beta_kj]
      new_expect_delta_beta <- expect_delta * mu_beta;
      model$params$expect_delta_beta[k, j] <- new_expect_delta_beta;
      model$params$expect_beta[k, j] <- new_expect_delta_beta;
      
      # E[(beta_kj)^2]
      model$params$expect_beta2[k, j] <- (1 - expect_delta) / tau + expect_delta * (mu_beta^2 + sigma2_beta);
      
      # Only target column j changes after this coordinate update
      model$params$expect_log_m[, j] <- expect_n + C_k * (new_expect_delta_beta - old_expect_delta_beta);
    }
  }
  
  # Refresh all m moments
  model <- update_m_moments(model);
  
  model
}

# Stirling lower-bound
z_shape_normalizer_lower_bound <- function(model) {
  with(model$params,
       -gamma * expect_m_log_m + 0.5 * expect_log_m + gamma * expect_m + 0.5 * log(gamma / (2 * pi)) - expect_m_inverse / (12 * gamma))
}

# Update q(z_ij, v_ij)
update_q_z_v <- function(model) {
  x <- model$data$x;
  u <- model$data$u;
  N <- model$data$N;
  J <- model$data$J;
  gamma <- model$params$gamma;
  rho <- model$params$rho;
  xi_z <- model$params$xi_z;
  
  # a_{z_{ij}}
  a_z <- x + gamma * model$params$expect_m;
  
  # b_{z_{ij}}
  b_z <- gamma + ((rho + x) * (1 - u)) / xi_z;
  
  if (
    any(!is.finite(a_z)) || any(a_z <= 0) ||
    any(!is.finite(b_z)) || any(b_z <= 0)
  ) {
    stop("Invalid Gamma parameters in q(z|v=1) update.");
  }
  
  # Store new q(z|v=1)
  model$vparams$z$a <- a_z;
  model$vparams$z$b <- b_z;
  
  # E[log(eta_j)]
  expect_log_eta_matrix <- matrix(model$params$expect_log_eta, nrow = N, ncol = J, byrow = TRUE);
  
  # v_ij = 1
  log_c_one <- expect_log_eta_matrix + z_shape_normalizer_lower_bound(model) - (rho + x) * log(xi_z) - ((rho + x) * (rho - xi_z)) / xi_z + lgamma(a_z) - a_z * log(b_z);
  
  expect_v <- matrix(1, N, J);
  zero_index <- which(x == 0, arr.ind=TRUE);
  
  if (nrow(zero_index) > 0) {
    for (r in seq_len(nrow(zero_index))) {
      i <- zero_index[r, 1];
      j <- zero_index[r, 2];
      
      # v_ij = 0
      log_c_zero <- model$params$expect_log_one_minus_eta[j] - rho * log(rho);
      
      # Normalization
      log_denominator <- log_sum_exp(c(log_c_zero, log_c_one[i, j]));
      
      # E[v_ij]
      expect_v[i, j] <- exp(log_c_one[i, j] - log_denominator);
    }
  }
  
  model$params$expect_v <- expect_v;
  
  model
}

# Update q(eta_j)
update_q_eta <- function(model) {
  model$vparams$eta <- list(
    a = model$params$alpha_eta + colSums(model$params$expect_v),
    b = model$params$beta_eta + model$data$N - colSums(model$params$expect_v)
  );
  
  model$params$expect_log_eta <- expect_log_x_beta(model$vparams$eta);
  model$params$expect_log_one_minus_eta <- expect_log_1mx_beta(model$vparams$eta);
  
  model
}

# Update q(psi_k)
update_q_psi <- function(model) {
  expect_delta <- model$params$expect_delta;
  xi_delta <- model$params$xi_delta;
  
  model$vparams$psi <- list(
    a = model$params$alpha_psi +
      rowSums(expect_delta + (1 - expect_delta) * (1 - xi_delta)),
    b = model$params$beta_psi +
      rowSums((1 - expect_delta) * xi_delta)
  );
  
  model$params$expect_log_psi <- expect_log_x_beta(model$vparams$psi);
  model$params$expect_log_one_minus_psi <- expect_log_1mx_beta(model$vparams$psi);
  
  model
}

# Update q(phi_j)
update_q_phi <- function(model) {
  expect_delta <- model$params$expect_delta;
  xi_delta <- model$params$xi_delta;
  
  model$vparams$phi <- list(
    a = model$params$alpha_phi + colSums(expect_delta),
    b = model$params$beta_phi +
      colSums((1 - expect_delta) * (1 - xi_delta))
  );
  
  model$params$expect_log_phi <- expect_log_x_beta(model$vparams$phi);
  model$params$expect_log_one_minus_phi <- expect_log_1mx_beta(model$vparams$phi);
  
  model
}

# Complete ELBO lower bound
elbo_components <- function(model, include_fixed_a_prior=TRUE) {
  x <- model$data$x;
  u <- model$data$u;
  N <- model$data$N;
  J <- model$data$J;
  K <- model$data$K;
  
  with(model$params, {
    a_z <- model$vparams$z$a;
    b_z <- model$vparams$z$b;
    expect_log_z_given_v_one <- digamma(a_z) - log(b_z);
    expect_z_given_v_one <- a_z / b_z;
    
    # E[log p(x|z)] + E[log p(z|m,v)] - E[log q(z|v)]
    observation_constant <- lgamma(x + rho) - lgamma(x + 1) - lgamma(rho) + rho * log(rho) + x_log_y(x, 1 - u);
    observation_auxiliary <- -(rho + x) * log(xi_z) - ((rho + x) / xi_z) * (rho + (1 - u) * expect_z_given_v_one - xi_z);
    latent_shape_lower_bound <- z_shape_normalizer_lower_bound(model);
    z_variational_part <- (x + gamma * expect_m - a_z) * expect_log_z_given_v_one - gamma * expect_z_given_v_one - a_z * log(b_z) + lgamma(a_z) + a_z;
    z_block <- sum(expect_v * (observation_constant + observation_auxiliary + latent_shape_lower_bound + z_variational_part));
    
    # E[log p(v|eta)] - E[log q(v)]
    v_eta_block <- sum(sweep(expect_v, 2, expect_log_eta, "*") + sweep(1 - expect_v, 2, expect_log_one_minus_eta, "*")) + entropy_bern_d(expect_v);
    
    # E[log p(beta)] - E[log q(beta|delta)]
    mu_beta <- model$vparams$beta_delta$mu_beta;
    sigma2_beta <- model$vparams$beta_delta$sigma2_beta;
    beta_conditional_block <- sum(-0.5 * expect_delta * (tau * (mu_beta^2 + sigma2_beta) - 1 - log(tau * sigma2_beta)));
    
    # -E[log q(delta)]
    delta_variational_block <- entropy_bern_d(expect_delta);
    
    # Jensen lower bound for E[log p(delta|psi,phi)]
    expect_log_psi_matrix <- matrix(expect_log_psi, nrow = K, ncol = J);
    expect_log_one_minus_psi_matrix <- matrix(expect_log_one_minus_psi, nrow = K, ncol = J);
    expect_log_phi_matrix <- matrix(expect_log_phi, nrow = K, ncol = J, byrow = TRUE);
    expect_log_one_minus_phi_matrix <- matrix(expect_log_one_minus_phi, nrow = K, ncol = J, byrow = TRUE);
    H_xi_delta <- bern_h(xi_delta);
    
    delta_prior_block <- sum(expect_delta * (expect_log_psi_matrix + expect_log_phi_matrix) +
                               (1 - expect_delta) * (xi_delta * expect_log_one_minus_psi_matrix +
                                                       (1 - xi_delta) * expect_log_psi_matrix +
                                                       (1 - xi_delta) * expect_log_one_minus_phi_matrix + H_xi_delta));
    
    beta_delta_block <- beta_conditional_block + delta_variational_block + delta_prior_block;
    
    # E[log p] - E[log q] for the Beta variational factors
    psi_block <- beta_prior_minus_q(model$vparams$psi, alpha_psi, beta_psi);
    phi_block <- beta_prior_minus_q(model$vparams$phi, alpha_phi, beta_phi);
    eta_block <- beta_prior_minus_q(model$vparams$eta, alpha_eta, beta_eta);
    
    # a is fixed in the first implementation, so there is no q(a) term
    a_block <- if (include_fixed_a_prior) {
      log_prior_a(model$data$a, alpha)
    } else {
      0
    };
    
    total <- z_block + v_eta_block + beta_delta_block + psi_block + phi_block + eta_block + a_block;
    
    list(
      z = z_block,
      v_eta = v_eta_block,
      beta_delta = beta_delta_block,
      psi = psi_block,
      phi = phi_block,
      eta = eta_block,
      a_fixed = a_block,
      total = total
    )
  })
}

# Variational evidence lower bound
elbo <- function(model) {
  elbo_components(model, include_fixed_a_prior=TRUE)$total
}

# Store the ELBO after one coordinate update
append_elbo_trace <- function(model, iteration, step, previous_elbo, allow_decrease=FALSE, verbose=1) {
  components <- elbo_components(model, include_fixed_a_prior=TRUE);
  current_elbo <- components$total;
  delta_elbo <- current_elbo - previous_elbo;
  
  trace_row <- data.frame(
    iteration = iteration,
    step = step,
    elbo = current_elbo,
    delta = delta_elbo,
    z = components$z,
    v_eta = components$v_eta,
    beta_delta = components$beta_delta,
    psi = components$psi,
    phi = components$phi,
    eta = components$eta,
    a_fixed = components$a_fixed,
    stringsAsFactors = FALSE
  );
  
  model$elbo_trace <- rbind(model$elbo_trace, trace_row);
  
  if (verbose >= 4) {
    if (!allow_decrease && is.finite(delta_elbo) && delta_elbo < 0) {
      message("WARN: ELBO decreased after ", step, ": ", delta_elbo);
    }
    
    if (allow_decrease && is.finite(delta_elbo) && delta_elbo < 0) {
      message("NOTE: ELBO decreased after ", step, " (Taylor q(beta,delta) update): ", delta_elbo
      );
    }
    
    message("delta elbo after ", step, ": ", delta_elbo);
  }
  
  list(model = model, elbo = current_elbo)
}