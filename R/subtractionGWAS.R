#' Vectorised SNP effects for square first-order factor models
#'
#' Estimates effects on k factors from k observed traits, using the fixed
#' loading and sampling covariance conventions of analytic userGWAS().
#' @param covstruc LDSC result containing S, V and I. S names at least two traits.
#' @param SNPs Data frame from sumstats(), with beta.<trait> and se.<trait>
#'   columns. Other columns are retained as metadata in their original order.
#' @param model First-order lavaan syntax with k factors and k observed traits,
#'   and one unconstrained '~ SNP' path per factor. X must be invertible.
#' @param usermod Optional no-SNP usermodel() result or its results table.
#'   The caller must establish convergence for supplied parameter tables.
#' @param batch_size Positive integer limiting the number of SNPs per batch.
#' @return A data frame with metadata, beta_<factor>, SE_<factor>,
#'   Z_beta_<factor>, p_val_<factor>, and subtraction_status. Invalid input rows
#'   retain their position with NA estimates. Q is omitted because the SNP
#'   model has zero residual degrees of freedom. The subtraction_diagnostics
#'   attribute records trait/factor order, loadings, adjusted intercepts,
#'   solver, condition number, invalid rows and no-SNP smoothing indicators.
#' @details
#' The no-SNP fit uses usermodel(), DWLS, std.lv=FALSE, CFIcalc=FALSE,
#' imp_cov=FALSE and fix_resid=TRUE, as in the analytic reference. Its S/V
#' smoothing is delegated to usermodel(); diagnostics flag non-positive
#' definite inputs. For supplied fits, smoothing history is unknown.
#' S and V must follow the same trait/vech ordering as an LDSC result.
#'
#' Each SNP covariance is diag(SE) I diag(SE), after raising univariate
#' intercepts below one to one. Named I axes are aligned to S. No SNP-variance
#' scaling, allele harmonisation, sign conversion or study-specific SE filter
#' is added. Only finite betas and finite positive SEs are estimated.
#'
#' The inverse transformation is prepared once using triangular solution for
#' lower-triangular X and a general linear solve otherwise. SNP arrays use
#' standard double precision. The existing condition policy warns at 1e6 and
#' rejects above 1e8. Non-finite results or non-positive variances stop the call.
#' Loadings remain fixed, as in analytic userGWAS(); their uncertainty is not
#' propagated. A square SNP model need not imply an exact no-SNP fit to S.
#'
#' Cholesky models use orthogonal unit-variance factors and zero observed
#' residual variances. Trait ordering defines the sequence of residuals.
#' General invertible square loadings are also supported. Latent regressions,
#' higher-order factors, multiple groups, constraints, defined parameters and
#' constrained SNP effects are unsupported. Use userGWAS() for other models.
#' This is a development implementation; userGWAS() remains unchanged.
#' @seealso userGWAS, usermodel, sumstats, ldsc
#' @export
subtractionGWAS <- function(covstruc, SNPs, model, usermod = NULL,
                           batch_size = 100000) {
  if (!is.numeric(batch_size) || length(batch_size) != 1L ||
      !is.finite(batch_size) || batch_size < 1 || batch_size != floor(batch_size))
    stop("batch_size must be a positive integer", call. = FALSE)
  prepared <- .subtraction_prepare_model(covstruc, model, usermod)
  traits <- prepared$traits
  factors <- prepared$factors
  SNPs <- as.data.frame(SNPs)
  if (anyDuplicated(names(SNPs))) stop("Duplicate SNP column names", call. = FALSE)
  required <- c(paste0("beta.", traits), paste0("se.", traits))
  if (!all(required %in% names(SNPs)))
    stop("SNPs must contain beta.<trait> and se.<trait> for every named trait",
         call. = FALSE)
  if (!all(vapply(SNPs[required], is.numeric, logical(1))))
    stop("GWAS beta and SE columns must be numeric", call. = FALSE)
  extra <- setdiff(grep("^(beta|se)\\.", names(SNPs), value = TRUE), required)
  if (length(extra)) stop("Unexpected beta/SE columns for other traits", call. = FALSE)
  output <- SNPs[setdiff(names(SNPs), required)]
  fields <- unlist(lapply(factors, function(f)
    paste0(c("beta_", "SE_", "Z_beta_", "p_val_"), f)), use.names = FALSE)
  reserved <- c(fields, "Q_omnibus", "Q_omnibus_df", "Q_omnibus_pval",
                "subtraction_status")
  if (any(reserved %in% names(output)))
    stop("Metadata column conflicts with a result column", call. = FALSE)
  for (field in fields) output[[field]] <- rep(NA_real_, nrow(SNPs))
  beta <- as.matrix(SNPs[paste0("beta.", traits)])
  se <- as.matrix(SNPs[paste0("se.", traits)])
  bad_beta <- rowSums(!is.finite(beta)) > 0
  bad_se <- rowSums(!is.finite(se) | se <= 0) > 0
  output$subtraction_status <- ifelse(bad_beta & bad_se, "nonfinite_beta;invalid_se",
    ifelse(bad_beta, "nonfinite_beta", ifelse(bad_se, "invalid_se", "ok")))
  valid <- which(!bad_beta & !bad_se)
  engine <- .subtraction_prepare_engine(prepared$loading, prepared$sampling)
  if (length(valid)) {
    for (start in seq(1, length(valid), by = batch_size)) {
      rows <- valid[seq.int(start, min(start + batch_size - 1, length(valid)))]
      result <- .subtraction_apply_engine(beta[rows, , drop = FALSE],
        se[rows, , drop = FALSE], engine)
      for (j in seq_along(factors)) {
        output[[paste0("beta_", factors[j])]][rows] <- result$beta[, j]
        output[[paste0("SE_", factors[j])]][rows] <- result$se[, j]
        output[[paste0("Z_beta_", factors[j])]][rows] <- result$z[, j]
        output[[paste0("p_val_", factors[j])]][rows] <- result$p[, j]
      }
    }
  }
  attr(output, "subtraction_diagnostics") <- list(
    method = "SubtractionSEM", implementation = "0.2.0",
    traits = traits, factors = factors, solver = engine$solver,
    smoothing = prepared$smoothing, loading = prepared$loading,
    sampling = prepared$sampling, condition_number = engine$condition_number,
    condition_warning = 1e6, condition_rejection = 1e8,
    invalid_rows = sum(bad_beta | bad_se), model_source = prepared$model_source,
    q_status = "not_applicable_just_identified_model")
  output
}
