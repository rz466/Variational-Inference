# Fit the ExprSim variational model
exprsim_vi <- function(data, hparams, model=NULL, tol=1e-2, max.iter=50, verbose=1) {
  
  if (is.null(model)) {
    fit <- exprsim_vi_local(
      data=data,
      hparams=hparams,
      verbose=verbose,
      tol=tol,
      max.iter=max.iter
    );
  } else {
    stopifnot(is(model, "exprsim_vi_fit"));
    fit <- exprsim_vi_local(
      model=model,
      verbose=verbose,
      tol=tol,
      max.iter=max.iter
    );
  }
  
  fit
}

# Initialize an exprsim_vi_fit object
make_model <- function(data, params0=NULL, hparams) {
  
  if (!is(data, "exprsim_vi_data")) {
    data <- prepare_data(data);
  }
  
  prepared_params <- prepare_params(params0, data, hparams);
  
  # Combine params and hparams
  model <- list(
    data = data,
    params = c(prepared_params$params, hparams),
    vparams = prepared_params$vparams,
    elbo = numeric(0),
    elbo_trace = data.frame()
  );
  
  # Expected log Beta quantities
  model$params$expect_log_psi <- expect_log_x_beta(model$vparams$psi);
  model$params$expect_log_one_minus_psi <- expect_log_1mx_beta(model$vparams$psi);
  
  model$params$expect_log_phi <- expect_log_x_beta(model$vparams$phi);
  model$params$expect_log_one_minus_phi <- expect_log_1mx_beta(model$vparams$phi);
  
  model$params$expect_log_eta <- expect_log_x_beta(model$vparams$eta);
  model$params$expect_log_one_minus_eta <- expect_log_1mx_beta(model$vparams$eta);
  
  # Initialize auxiliary quantities and variational blocks
  model <- update_xi_delta(model);
  model <- update_m_moments(model);
  model <- update_q_z_v(model);
  model <- update_xi_z(model);
  model <- update_q_eta(model);
  
  structure(model, class="exprsim_vi_fit")
}

# Predict expression counts
# Predict expression counts
predict.exprsim_vi_fit <- function(object, u_new, a_new) {
  
  N <- nrow(a_new);
  J <- object$data$J;
  K <- object$data$K;
  
  data_new <- list(u = u_new, a = a_new, N = N, J = J, K = K, C = (1 - u_new[, seq_len(K), drop=FALSE]) * a_new);
  
  m_moments <- compute_m_moments(data=data_new, expect_delta=object$params$expect_delta,
                                 mu_beta=object$vparams$beta_delta$mu_beta, sigma2_beta=object$vparams$beta_delta$sigma2_beta);
  
  expect_eta <- expect_x_beta(object$vparams$eta);
  
  expect_x <- (1 - u_new) * m_moments$expect_m * matrix(expect_eta, nrow=N, ncol=J, byrow=TRUE);
  
  expect_x
}

# Extract edge probabilities or beta means
coef.exprsim_vi_fit <- function(object, par=c("delta", "beta_mu", "beta_mean")) {
  par <- match.arg(par);
  
  if (par == "delta") {
    return(object$params$expect_delta);
  }
  
  if (par == "beta_mu") {
    return(object$vparams$beta_delta$mu_beta);
  }
  
  object$params$expect_beta
}

# Print fitted model
print.exprsim_vi_fit <- function(object) {
  niters <- length(object$elbo);
  
  elbo_value <- if (niters > 0) {
    object$elbo[niters]
  } else {
    elbo(object)
  }
  
  n_called <- sum(object$params$expect_delta > 0.5);
  
  cat("\nExprSim variational model\n\n");
  cat("ELBO lower bound: ", elbo_value, ",  niters: ", niters, "\n", sep="");
  cat("edges: ", n_called, " called out of ", object$data$K * object$data$J, "\n", sep="");
}