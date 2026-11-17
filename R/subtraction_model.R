# Internal model preparation. Uses the public usermodel() API; does not change
# the preparation or numerical behavior of the existing userGWAS() branches.
.subtraction_named_matrix <- function(x, traits, name) {
  if (!is.matrix(x) || !is.numeric(x) || !identical(dim(x), c(length(traits), length(traits))) ||
      any(!is.finite(x))) stop(name, " must be a finite numeric k by k matrix",
                             call. = FALSE)
  dn <- dimnames(x)
  if (!is.null(dn)) {
    if (length(dn[[1]]) != length(traits) || length(dn[[2]]) != length(traits) ||
        anyDuplicated(dn[[1]]) || anyDuplicated(dn[[2]]) ||
        !setequal(dn[[1]], traits) || !setequal(dn[[2]], traits))
      stop(name, " axis names must match S traits", call. = FALSE)
    x <- x[traits, traits, drop = FALSE]
  }
  if (!isTRUE(all.equal(x, t(x), tolerance = 1e-12, check.attributes = FALSE)))
    stop(name, " must be symmetric", call. = FALSE)
  # Resolve sub-tolerance asymmetry deterministically before the kernel uses I12.
  x <- (x + t(x)) / 2
  dimnames(x) <- list(traits, traits)
  x
}

.subtraction_prepare_model <- function(covstruc, model, usermod) {
  if (!is.list(covstruc) || is.null(covstruc$S) || is.null(covstruc$I))
    stop("covstruc must contain S and I", call. = FALSE)
  traits <- colnames(covstruc$S)
  if (length(traits) < 2 || anyNA(traits) || any(!nzchar(traits)) ||
      anyDuplicated(traits) || "SNP" %in% traits)
    stop("S must name at least two distinct observed traits", call. = FALSE)
  S <- covstruc$S
  if (is.null(rownames(S))) rownames(S) <- traits
  S <- .subtraction_named_matrix(S, traits, "S")
  if (any(diag(S) <= 0)) stop("S must have positive variances", call. = FALSE)
  sampling <- .subtraction_named_matrix(covstruc$I, traits, "I")
  diag(sampling) <- pmax(diag(sampling), 1)
  if (min(eigen(sampling, symmetric = TRUE, only.values = TRUE)$values) <= 0)
    stop("Adjusted I must be positive definite", call. = FALSE)
  if (!is.character(model) || length(model) != 1L || is.na(model) || !nzchar(model))
    stop("model must be one nonempty lavaan syntax string", call. = FALSE)
  # Include sem()'s automatic variances/covariances; auto=FALSE would turn
  # omitted free variances into fixed zeros when reconstructing the model.
  pt <- lavaan::lavaanify(model, auto = TRUE, auto.fix.first = TRUE)
  if (any(!pt$op %in% c("=~", "~~", "~")))
    stop("Constraints, intercepts and defined parameters are not supported",
         call. = FALSE)
  if (any(pt$group != 1L) || any(pt$block != 1L))
    stop("Multiple groups or blocks are unsupported", call. = FALSE)
  loadings <- pt[pt$op == "=~", , drop = FALSE]
  factors <- unique(loadings$lhs)
  if (length(factors) != length(traits) || any(factors %in% c(traits, "SNP")) ||
      !setequal(unique(loadings$rhs), traits))
    stop("Expected k first-order factors loading only on the k S traits",
         call. = FALSE)
  if (any(!pt$lhs %in% c(traits, factors, "SNP")) ||
      any(!pt$rhs %in% c(traits, factors, "SNP")))
    stop("Model contains an unknown variable", call. = FALSE)
  paths <- pt[pt$op == "~", , drop = FALSE]
  if (nrow(paths) != length(factors) || !setequal(paths$lhs, factors) ||
      any(paths$rhs != "SNP") || any(paths$free == 0) ||
      anyDuplicated(paths$free))
    stop("Specify one free '~ SNP' path for each factor; other regressions are unsupported",
         call. = FALSE)
  for (label in paths$label[nzchar(paths$label)])
    if (sum(pt$label == label) != 1L)
      stop("SNP effects must not share parameter labels", call. = FALSE)
  snp_cov <- pt[pt$user == 1 & pt$op == "~~" &
                  (pt$lhs == "SNP" | pt$rhs == "SNP"), ]
  if (nrow(snp_cov) && any(snp_cov$lhs != "SNP" | snp_cov$rhs != "SNP" |
                          snp_cov$free != 0 | snp_cov$ustart != 0))
    stop("Additional SNP covariance parameters are unsupported", call. = FALSE)
  # No measurement error may covary with a latent factor in this first version.
  cross <- pt$op == "~~" & ((pt$lhs %in% factors & pt$rhs %in% traits) |
                            (pt$rhs %in% factors & pt$lhs %in% traits))
  if (any(cross)) stop("Observed-latent residual covariances are unsupported",
                       call. = FALSE)

  model_source <- if (is.null(usermod)) "usermodel_DWLS" else "supplied_usermod"
  smoothing <- list(S = NA, V = NA, source = "supplied_fit_not_reconstructed")
  if (is.null(usermod)) {
    V <- covstruc$V
    if (!is.matrix(V) || !is.numeric(V) || !identical(dim(V), rep(as.integer(length(traits) * (length(traits) + 1L) / 2L), 2)) ||
        any(!is.finite(V))) stop("Fitting requires a finite k(k+1)/2 square V", call. = FALSE)
    if (!isTRUE(all.equal(V, t(V), tolerance = 1e-12, check.attributes = FALSE)))
      stop("V must be symmetric", call. = FALSE)
    smoothing <- list(S = min(eigen(S, symmetric = TRUE, only.values = TRUE)$values) <= 0,
      V = min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) <= 0,
      source = "usermodel_nearPD")
    covstruc$S <- S
    # Reconstruct measurement syntax from the parsed table so semicolons and
    # '~ SNP' on a shared line cannot accidentally remove measurement terms.
    measurement <- pt[pt$lhs != "SNP" & pt$rhs != "SNP", , drop = FALSE]
    terms <- vapply(seq_len(nrow(measurement)), function(i) {
      row <- measurement[i, ]
      modifier <- if (row$free == 0) paste0(format(row$ustart, digits = 17), "*")
                  else if (nzchar(row$label)) paste0(row$label, "*") else "NA*"
      paste(row$lhs, row$op, paste0(modifier, row$rhs))
    }, character(1))
    # usermodel() reads V and S positionally; use a canonical ordering here.
    fit_input <- list(V = V, S = S, I = sampling)
    usermod <- tryCatch(usermodel(fit_input, estimation = "DWLS",
      model = paste(terms, collapse = "\n"), CFIcalc = FALSE,
      std.lv = FALSE, imp_cov = FALSE, fix_resid = TRUE),
      error = function(e) stop("No-SNP model fitting failed: ", conditionMessage(e),
                               call. = FALSE))
    if (is.null(usermod)) stop("The no-SNP model did not return a valid fit",
                              call. = FALSE)
  }
  table <- if (is.data.frame(usermod)) usermod else usermod$results
  if (!is.data.frame(table) || !all(c("lhs", "op", "rhs") %in% names(table)))
    stop("usermod must contain a no-SNP parameter table", call. = FALSE)
  if (!"Unstand_Est" %in% names(table) && "Unstandardized_Estimate" %in% names(table))
    table$Unstand_Est <- table$Unstandardized_Estimate
  if (!is.numeric(table$Unstand_Est) || any(!is.finite(table$Unstand_Est)))
    stop("Parameter estimates must be finite numeric values", call. = FALSE)
  keys <- paste(table$lhs, table$op, table$rhs)
  if (anyDuplicated(keys)) stop("Duplicate parameter rows in usermod", call. = FALSE)
  if (any(!table$op %in% c("=~", "~~")) ||
      any(!table$lhs %in% c(traits, factors)) ||
      any(!table$rhs %in% c(traits, factors)))
    stop("usermod must describe this no-SNP first-order model", call. = FALSE)
  if (any(table$op == "~~" & table$lhs == table$rhs & table$Unstand_Est < 0))
    stop("Negative model variance", call. = FALSE)
  # Fixed rows can be absent in a supplied free-parameter table; fill them
  # from parsed syntax, while requiring every free measurement parameter.
  measurement <- pt[pt$lhs != "SNP" & pt$rhs != "SNP", , drop = FALSE]
  for (i in seq_len(nrow(measurement))) {
    row <- measurement[i, ]
    match_row <- table$op == row$op & table$lhs == row$lhs & table$rhs == row$rhs
    if (row$op == "~~") match_row <- match_row |
      (table$op == "~~" & table$lhs == row$rhs & table$rhs == row$lhs)
    hit <- table[match_row, ]
    if (nrow(hit) > 1L) stop("Duplicate symmetric covariance rows", call. = FALSE)
    if (!nrow(hit)) {
      if (row$free != 0) stop("Missing estimated loading or covariance: ",
        row$lhs, " ", row$op, " ", row$rhs, call. = FALSE)
      added <- table[NA_integer_, , drop = FALSE]
      added$lhs <- row$lhs; added$op <- row$op; added$rhs <- row$rhs
      added$Unstand_Est <- row$ustart
      table <- rbind(table, added)
    } else if (row$free == 0 &&
               !isTRUE(all.equal(hit$Unstand_Est, row$ustart, tolerance = 1e-12))) {
      stop("Supplied parameter disagrees with a fixed model parameter", call. = FALSE)
    }
  }
  variances <- table[table$op == "~~" & table$lhs == table$rhs, ]
  if (any(variances$Unstand_Est < 0)) stop("Negative model variance", call. = FALSE)
  covariances <- table[table$op == "~~" & table$lhs != table$rhs, ]
  for (i in seq_len(nrow(covariances))) {
    row <- covariances[i, ]
    v <- variances$Unstand_Est[match(c(row$lhs, row$rhs), variances$lhs)]
    if (anyNA(v) || abs(row$Unstand_Est) > sqrt(prod(v)) + 1e-12)
      stop("Missing variances or inadmissible model covariance", call. = FALSE)
  }
  x <- matrix(0, length(traits), length(factors), dimnames = list(traits, factors))
  fitted_loadings <- table[table$op == "=~", ]
  declared <- paste(loadings$lhs, loadings$rhs)
  if (any(!paste(fitted_loadings$lhs, fitted_loadings$rhs) %in% declared))
    stop("usermod loadings do not match the declared model", call. = FALSE)
  for (i in seq_len(nrow(loadings))) {
    row <- loadings[i, ]
    hit <- table[table$lhs == row$lhs & table$rhs == row$rhs & table$op == "=~", ]
    if (!nrow(hit) && row$free != 0)
      stop("Missing estimated loading: ", row$lhs, " =~ ", row$rhs, call. = FALSE)
    value <- if (nrow(hit)) hit$Unstand_Est else row$ustart
    if (row$free == 0 && !isTRUE(all.equal(value, row$ustart, tolerance = 1e-12)))
      stop("Supplied loading disagrees with a fixed model loading", call. = FALSE)
    x[row$rhs, row$lhs] <- value
  }
  list(loading = x, sampling = sampling, traits = traits, factors = factors,
       model_source = model_source, smoothing = smoothing)
}
