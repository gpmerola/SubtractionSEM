# Canonical general engine shared by the package and research runner.
# All SNP work is vectorised; matrices depending only on the model are prepared once.
.subtraction_prepare_engine <- function(loading, sampling) {
  k <- nrow(loading)
  condition <- kappa(loading, exact = TRUE)
  if (!is.finite(condition)) stop("Singular loading matrix", call. = FALSE)
  if (condition > 1e8)
    stop("Loading-matrix condition number exceeds rejection threshold 1e8", call. = FALSE)
  if (condition >= 1e6)
    warning("Loading-matrix condition number exceeds 1e6", call. = FALSE)
  triangular <- all(loading[upper.tri(loading)] == 0)
  inverse <- if (triangular) forwardsolve(loading, diag(k)) else solve(loading)
  coefficients <- lapply(seq_len(k), function(f) tcrossprod(inverse[f, ]) * sampling)
  list(inverse = inverse, coefficients = coefficients, condition_number = condition,
       solver = if (triangular) "triangular" else "general")
}

.subtraction_apply_engine <- function(beta, se, engine) {
  n <- nrow(beta); k <- ncol(beta)
  estimates <- errors <- matrix(0, n, k)
  # Extract columns once per batch; repeated matrix slices allocate fresh vectors.
  beta_columns <- lapply(seq_len(k), function(a) beta[, a])
  se_columns <- lapply(seq_len(k), function(a) se[, a])
  se_squared <- lapply(se_columns, function(x) x^2)
  for (f in seq_len(k)) {
    estimate <- variance <- numeric(n)
    coefficients <- engine$coefficients[[f]]
    for (a in seq_len(k)) {
      # Exact zeros only: no thresholding or approximation of small coefficients.
      if (engine$inverse[f, a] != 0)
        estimate <- estimate + engine$inverse[f, a] * beta_columns[[a]]
      if (coefficients[a, a] != 0)
        variance <- variance + coefficients[a, a] * se_squared[[a]]
      if (a < k) for (b in seq.int(a + 1L, k)) {
        if (coefficients[a, b] != 0)
          variance <- variance + 2 * coefficients[a, b] * se_columns[[a]] * se_columns[[b]]
      }
    }
    if (any(!is.finite(variance) | variance <= 0))
      stop("Non-positive or non-finite transformed variance", call. = FALSE)
    estimates[, f] <- estimate
    errors[, f] <- sqrt(variance)
  }
  z <- estimates / errors
  if (any(!is.finite(estimates)) || any(!is.finite(z)))
    stop("Non-finite transformed result", call. = FALSE)
  list(beta = estimates, se = errors, z = z, p = 2 * pnorm(-abs(z)))
}
