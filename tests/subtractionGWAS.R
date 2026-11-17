# Synthetic-only integration tests, run by R CMD check or Rscript.
library(GenomicSEM)
library(SubtractionSEM)
set.seed(9217)
checks <- 0L
check <- function(value) {
  stopifnot(isTRUE(value))
  checks <<- checks + 1L
}
reject <- function(expr, pattern) {
  error <- tryCatch({force(expr); NULL}, error = function(e) conditionMessage(e))
  if (is.null(error) || !grepl(pattern, error))
    stop("Expected ", pattern, "; received ", if (is.null(error)) "success" else error,
         " in ", paste(deparse(substitute(expr)), collapse = " "))
  check(!is.null(error) && grepl(pattern, error))
}
near <- function(a, b, atol = 5e-13, rtol = 5e-12) {
  check(identical(dim(a), dim(b)) && length(a) == length(b) &&
          all(is.finite(a)) && all(is.finite(b)) &&
          all(abs(a - b) <= atol + rtol * abs(b)))
}
traits <- c("T1", "T2")
loading <- matrix(c(.8, 0, .6, .7), 2, byrow = TRUE,
                  dimnames = list(traits, c("Common", "Unique")))
S <- tcrossprod(loading)
I <- matrix(c(1.1, .2, .2, 1), 2, dimnames = list(traits, traits))
covstruc <- list(V = diag(.001, 3), S = S, I = I)
model <- paste("Common =~ NA*T1 + T2", "Unique =~ NA*T2",
  "Common ~~ 1*Common", "Unique ~~ 1*Unique", "Common ~~ 0*Unique",
  "T1 ~~ 0*T1", "T2 ~~ 0*T2", "Common ~ SNP", "Unique ~ SNP", sep = "\n")
params <- data.frame(lhs = c("Common", "Common", "Unique", "Common", "Unique"),
  op = c("=~", "=~", "=~", "~~", "~~"),
  rhs = c("T1", "T2", "T2", "Common", "Unique"),
  Unstand_Est = c(.8, .6, .7, 1, 1))
snps <- data.frame(SNP = paste0("synthetic_", seq_len(40)), CHR = 1L,
  BP = seq_len(40), MAF = .2, A1 = "A", A2 = "C",
  beta.T1 = rnorm(40, sd = .03), beta.T2 = rnorm(40, sd = .02),
  se.T1 = runif(40, .01, .04), se.T2 = runif(40, .01, .04))
estimate <- function(SNPs = snps, fit = params, syntax = model, cov = covstruc, batch = 7)
  subtractionGWAS(cov, SNPs, syntax, fit, batch_size = batch)
fields <- c("beta_Common", "SE_Common", "Z_beta_Common", "p_val_Common",
            "beta_Unique", "SE_Unique", "Z_beta_Unique", "p_val_Unique")
result <- estimate()
check(identical(result$SNP, snps$SNP))
check(all(result$subtraction_status == "ok"))
check(!any(c("Q_omnibus", "Q_omnibus_df", "Q_omnibus_pval") %in% names(result)))
check(identical(result[fields], estimate(batch = 1)[fields]))
check(identical(result[fields], estimate(batch = 1000)[fields]))
check(!any(grepl("subtraction_.*cpp$", getNamespaceExports("SubtractionSEM"))))

# Reference path, including its actual output column names and unchanged API.
official <- suppressWarnings(userGWAS(covstruc, snps, model = model,
                     analytic = TRUE, usermod = params, batch_size = 13))
for (field in fields) near(result[[field]], official[[field]],
                          atol = if (startsWith(field, "Z_")) 5e-12 else 5e-13)

# Relative/log-scale agreement also detects discrepancies hidden by tiny P values.
check_p <- function(a, b) {
  check(identical(is.na(a), is.na(b)))
  check(identical(a == 0, b == 0))
  positive <- is.finite(a) & is.finite(b) & a > 0 & b > 0
  check(all(abs(a[positive]-b[positive])/pmax(a[positive],b[positive]) < 1e-10))
  check(all(abs(log10(a[positive])-log10(b[positive])) < 1e-10))
}
for (f in c("Common", "Unique")) check_p(result[[paste0("p_val_",f)]], official[[paste0("p_val_",f)]])
tiny <- snps[1:3, ]
tiny$beta.T1 <- c(8, 20, 30) * tiny$se.T1
small <- estimate(SNPs=tiny)
reference_small <- suppressWarnings(userGWAS(covstruc, tiny, model=model, analytic=TRUE, usermod=params))
check(min(small$p_val_Common) < 1e-150)
check_p(small$p_val_Common, reference_small$p_val_Common)

# Independent matrix implementation across covariance scales and loading signs.
for (rho in c(-.75, 0, .75, .95)) {
  cov <- covstruc
  cov$I <- matrix(c(1, rho, rho, 1), 2)
  for (sign in c(-1, 1)) {
    fit <- params
    fit$Unstand_Est[3] <- sign * .7
    x <- loading
    x[2, 2] <- sign * .7
    out <- estimate(fit = fit, cov = cov)
    inv <- solve(x)
    beta <- t(inv %*% t(as.matrix(snps[c("beta.T1", "beta.T2")])))
    se <- t(vapply(seq_len(nrow(snps)), function(i) {
      d <- diag(c(snps$se.T1[i], snps$se.T2[i]))
      sqrt(diag(inv %*% d %*% cov$I %*% d %*% t(inv)))
    }, numeric(2)))
    for (j in seq_len(2)) {
      f <- colnames(loading)[j]
      near(out[[paste0("beta_", f)]], beta[, j])
      near(out[[paste0("SE_", f)]], se[, j])
      near(out[[paste0("Z_beta_", f)]], beta[, j] / se[, j], atol = 5e-12)
      near(out[[paste0("p_val_", f)]], 2 * pnorm(-abs(beta[, j] / se[, j])))
    }
  }
}
check(identical(result[fields], estimate(snps[rev(names(snps))])[fields]))
reordered <- covstruc
reordered$S <- S[2:1, 2:1]
for (field in fields) near(result[[field]], estimate(cov = reordered)[[field]], atol = 5e-12)
reordered$I <- I[2:1, 2:1]
for (field in fields) near(result[[field]], estimate(cov = reordered)[[field]], atol = 5e-12)
clamped <- covstruc
clamped$I[2, 2] <- .8
check(identical(result[fields], estimate(cov = clamped)[fields]))
small <- snps
small$se.T1[1] <- 1e-7
check(estimate(small)$subtraction_status[1] == "ok")
invalid <- snps
invalid$beta.T1[2] <- NA_real_
invalid$se.T2[3] <- 0
invalid$se.T1[4] <- -1
invalid$beta.T2[5] <- Inf
invalid$se.T2[5] <- Inf
out <- estimate(invalid)
check(identical(out$SNP, invalid$SNP))
check(all(is.na(as.matrix(out[2:5, fields]))))
check(identical(out$subtraction_status[2:5],
                c("nonfinite_beta", "invalid_se", "invalid_se", "nonfinite_beta;invalid_se")))
check(nrow(estimate(snps[FALSE, ])) == 0L)
check(nrow(estimate(snps[1, ])) == 1L)
check(all(is.na(as.matrix(estimate(invalid[2:5, ])[fields]))))
check(identical(estimate(fit = list(results = params))[fields], result[fields]))
# Omitted latent variances are free sem() parameters, not fixed zeros.
default_variance_model <- gsub("Common ~~ 1\\*Common\n|Unique ~~ 1\\*Unique\n", "", model)
check(identical(estimate(syntax = default_variance_model)[fields], result[fields]))
fixed_model <- sub("Unique =~ NA\\*T2", "Unique =~ 0.7*T2", model)
check(identical(estimate(fit = params[-3, ], syntax = fixed_model)[fields], result[fields]))
reject(estimate(syntax = sub("Unique =~ NA\\*T2", "Unique =~ 0.8*T2", model)),
       "fixed model parameter")

# Syntax and model errors must not silently select a different estimator.
reject(estimate(batch = 0), "positive integer")
reject(estimate(snps[setdiff(names(snps), "beta.T1")]), "must contain")
bad <- snps; bad$beta.T1 <- as.character(bad$beta.T1)
reject(estimate(bad), "numeric")
bad <- snps; bad$beta.T3 <- 0
reject(estimate(bad), "other traits")
bad <- snps; bad$subtraction_status <- "old"
reject(estimate(bad), "conflicts")
reject(estimate(syntax = paste(model, "Common ~ Unique", sep = "\n")), "regressions")
reject(estimate(syntax = sub("Common ~ SNP", "Common ~ 0*SNP", model)), "free")
reject(estimate(syntax = gsub("~ SNP", "~ a*SNP", model)), "Constraints|labels")
reject(estimate(syntax = paste(model, "New =~ Common", sep = "\n")), "first-order")
reject(estimate(syntax = paste(model, "x := 1", sep = "\n")), "defined parameters")
reject(estimate(fit = rbind(params, params[1, ])), "Duplicate")
reject(estimate(fit = params[-1, ]), "Missing estimated loading")
bad <- params; bad$Unstand_Est[4] <- -1
reject(estimate(fit = bad), "Negative")
bad <- params; bad$Unstand_Est[3] <- 0
reject(estimate(fit = bad), "Singular")
bad <- params; bad$Unstand_Est[3] <- 1e-10
reject(estimate(fit = bad), "rejection threshold")
bad <- covstruc; bad$I[1, 2] <- .8
reject(estimate(cov = bad), "symmetric")
bad <- covstruc; bad$I[,] <- matrix(c(1, 2, 2, 1), 2)
reject(estimate(cov = bad), "positive definite")
bad <- covstruc; dimnames(bad$I) <- list(c("A", "B"), c("A", "B"))
reject(estimate(cov = bad), "axis names")
bad <- params; bad$Unstand_Est[3] <- 1e-7
warned <- FALSE
withCallingHandlers(estimate(fit = bad), warning = function(w) {
  if (grepl("condition number", conditionMessage(w))) warned <<- TRUE
  invokeRestart("muffleWarning")
}) -> ignored
check(warned)

# Actual no-SNP fitting, not just a supplied parameter table.
fitted <- subtractionGWAS(covstruc, snps, model, batch_size = 11)
check(attr(fitted, "subtraction_diagnostics")$model_source == "usermodel_DWLS")
fit <- usermodel(covstruc, estimation = "DWLS",
  model = sub("\nCommon ~ SNP\nUnique ~ SNP", "", model),
  CFIcalc = FALSE, std.lv = FALSE, imp_cov = FALSE)
check(is.data.frame(fit$results))
supplied <- estimate(fit = fit)
for (field in fields) near(fitted[[field]], supplied[[field]], atol = 1e-9, rtol = 1e-8)
order_free <- subtractionGWAS(covstruc[c("S", "I", "V")], snps, model)
for (field in fields) near(fitted[[field]], order_free[[field]], atol = 1e-9, rtol = 1e-8)
semicolon <- gsub("\n", "; ", model)
for (field in fields) near(subtractionGWAS(covstruc, snps, semicolon)[[field]],
                           fitted[[field]], atol = 1e-9, rtol = 1e-8)
cat("SubtractionSEM integration checks passed:", checks, "\n")

# General square models: real public fits, not only kernel calls.
# Synthetic example generator; no empirical or participant data.
make_square_fixture <- function(k, dense = FALSE, n = 80L) {
  tr <- paste0("T", seq_len(k)); fs <- paste0("F", seq_len(k))
  x <- matrix(0, k, k, dimnames = list(tr, fs))
  x[lower.tri(x)] <- .15
  diag(x) <- seq(.6, .9, length.out = k)
  if (dense) x[upper.tri(x)] <- -.08
  cov <- list(V = diag(.0001, k*(k+1)/2), S = tcrossprod(x),
              I = matrix(.12, k, k, dimnames = list(tr, tr)))
  diag(cov$I) <- seq(.9, 1.2, length.out = k)
  syntax <- character(); rows <- list()
  for (f in seq_len(k)) {
    indices <- if (dense) seq_len(k) else seq.int(f, k)
    terms <- if (dense) paste0(x[indices, f], "*", tr[indices]) else
      c(paste0("NA*", tr[indices[1]]), tr[indices[-1]])
    syntax <- c(syntax, paste(fs[f], "=~", paste(terms, collapse=" + ")),
                paste(fs[f], "~~ 1*", fs[f]), paste(fs[f], "~ SNP"))
    for (i in indices) rows[[length(rows)+1L]] <- data.frame(
      lhs=fs[f], op="=~", rhs=tr[i], Unstand_Est=x[i,f])
    rows[[length(rows)+1L]] <- data.frame(lhs=fs[f],op="~~",rhs=fs[f],Unstand_Est=1)
    if (f < k) for (g in seq.int(f+1L,k))
      syntax <- c(syntax, paste(fs[f], "~~ 0*", fs[g]))
  }
  syntax <- c(syntax, paste(tr,"~~ 0*",tr))
  dat <- data.frame(SNP=paste0("general_",seq_len(n)),CHR=1L,BP=seq_len(n),
                    MAF=.25,A1="A",A2="G")
  for (i in seq_len(k)) dat[[paste0("beta.",tr[i])]] <- rnorm(n,sd=.025)
  for (i in seq_len(k)) dat[[paste0("se.",tr[i])]] <- runif(n,.01,.03)
  list(cov=cov, model=paste(syntax,collapse="\n"), params=do.call(rbind,rows),
       snps=dat, x=x, factors=fs)
}
for (k in c(2L,3L,4L,6L,12L)) {
  fixture <- make_square_fixture(k)
  output <- subtractionGWAS(fixture$cov,fixture$snps,fixture$model,fixture$params,batch_size=17)
  gold <- suppressWarnings(userGWAS(fixture$cov,fixture$snps,model=fixture$model,
                                    analytic=TRUE,usermod=fixture$params))
  for (f in fixture$factors) {
    for (prefix in c("beta_","SE_","Z_beta_","p_val_")) {
      field <- paste0(prefix,f)
      near(output[[field]],gold[[field]],atol=if(prefix=="Z_beta_") 5e-12 else 5e-13)
    }
    check_p(output[[paste0("p_val_",f)]],gold[[paste0("p_val_",f)]])
  }
  fitted <- subtractionGWAS(fixture$cov,fixture$snps,fixture$model)
  reference <- suppressWarnings(userGWAS(fixture$cov,fixture$snps,
                              model=fixture$model,analytic=TRUE))
  for (f in fixture$factors) for (prefix in c("beta_","SE_","Z_beta_","p_val_")) {
    field <- paste0(prefix,f)
    near(fitted[[field]],reference[[field]],atol=1e-8,rtol=1e-7)
  }
  reordered <- fixture$snps[rev(names(fixture$snps))]
  reordered <- subtractionGWAS(fixture$cov,reordered,fixture$model,fixture$params)
  for (f in fixture$factors) near(output[[paste0("beta_",f)]],reordered[[paste0("beta_",f)]])
  invalid <- fixture$snps; invalid[[paste0("se.T",k)]][2] <- NA_real_
  invalid[[paste0("beta.T",k)]][3] <- Inf
  bad <- subtractionGWAS(fixture$cov,invalid,fixture$model,fixture$params)
  check(all(bad$subtraction_status[2:3] != "ok"))
  check(all(is.na(bad$beta_F1[2:3])))
  check(nrow(subtractionGWAS(fixture$cov,invalid[FALSE,],fixture$model,fixture$params))==0)
  # Dense, invertible X is supported, not only triangular loadings.
  dense <- make_square_fixture(k,TRUE)
  out <- subtractionGWAS(dense$cov,dense$snps,dense$model,dense$params)
  expected <- as.matrix(dense$snps[paste0("beta.T",seq_len(k))]) %*% t(solve(dense$x))
  for (f in seq_len(k)) near(out[[paste0("beta_F",f)]],expected[,f])
}
cat("All square-model integration checks passed:",checks,"\n")

# Last Cholesky residual is invariant to permuting the preceding set.
f <- make_square_fixture(3)
a <- SubtractionSEM::subtractionGWAS(f$cov,f$snps,f$model,f$params)
order <- c(2L,1L,3L)
x <- t(chol(f$cov$S[order,order])); tr <- colnames(f$cov$S)[order]
y <- as.matrix(f$snps[paste0("beta.",tr)]) %*% t(solve(x))
near(a$beta_F3,y[,3])
# A model with fewer factors than traits must not silently change estimator.
reject(subtractionGWAS(f$cov,f$snps,"F =~ T1 + T2 + T3\nF ~ SNP",f$params),"first-order")
# General condition policy and singularity apply beyond 2x2.
reject(SubtractionSEM:::.subtraction_prepare_engine(diag(c(1,1,0)),diag(3)),"Singular")
reject(SubtractionSEM:::.subtraction_prepare_engine(diag(c(1,1,1e-10)),diag(3)),"rejection threshold")
# S/V repair is delegated to the same upstream fitter, explicitly diagnosed.
f <- make_square_fixture(3,TRUE)
f$cov$S[1,2] <- f$cov$S[2,1] <- .9
reject(suppressWarnings(subtractionGWAS(f$cov,f$snps,f$model)), "valid fit")
# A singular sampling V that upstream can repair remains supported.
f <- make_square_fixture(3)
f$cov$V[6,6] <- 0
smoothed <- suppressWarnings(subtractionGWAS(f$cov,f$snps,f$model))
check(isTRUE(attr(smoothed,"subtraction_diagnostics")$smoothing$V))
reference <- suppressWarnings(userGWAS(f$cov,f$snps,model=f$model,analytic=TRUE))
for (factor in f$factors) near(smoothed[[paste0("SE_",factor)]],reference[[paste0("SE_",factor)]],atol=1e-8,rtol=1e-7)
cat("General integration, ordering and smoothing checks passed:",checks,"\n")

check(is.na(attr(result,"subtraction_diagnostics")$smoothing$S))
asymmetric <- covstruc; asymmetric$V[1,2] <- .001
reject(subtractionGWAS(asymmetric,snps,model),"V must be symmetric")
cat("Final package integration checks passed:",checks,"\n")

# Regression oracle: preserve the original term order, including tiny nonzero
# coefficients. This frozen test implementation is not a runtime fallback.
legacy_kernel <- function(beta, se, engine) {
  n <- nrow(beta); k <- ncol(beta)
  estimates <- errors <- matrix(0, n, k)
  for (f in seq_len(k)) {
    variance <- numeric(n)
    for (a in seq_len(k)) {
      estimates[,f] <- estimates[,f] + engine$inverse[f,a] * beta[,a]
      variance <- variance + engine$coefficients[[f]][a,a] * se[,a]^2
      if (a < k) for (b in seq.int(a+1L,k))
        variance <- variance + 2 * engine$coefficients[[f]][a,b] * se[,a] * se[,b]
    }
    errors[,f] <- sqrt(variance)
  }
  z <- estimates/errors
  list(beta=estimates,se=errors,z=z,p=2*pnorm(-abs(z)))
}
prepare_engine <- getFromNamespace('.subtraction_prepare_engine','SubtractionSEM')
apply_engine <- getFromNamespace('.subtraction_apply_engine','SubtractionSEM')
set.seed(20260920)
for (k in c(2,3,4,6,12)) for (kind in c('triangular','dense','diagonal','tiny')) {
  loading <- matrix(rnorm(k*k,sd=.03),k,k) + diag(k)
  if (kind %in% c('triangular','tiny')) loading[upper.tri(loading)] <- 0
  if (kind == 'diagonal') loading <- diag(seq(.6,1.2,length.out=k))
  if (kind == 'tiny') loading[1,k] <- 1e-20
  sampling <- matrix(.15,k,k); diag(sampling) <- 1
  engine <- prepare_engine(loading,sampling)
  beta <- matrix(rnorm(97*k,sd=.05),97,k)
  se <- matrix(runif(97*k,.01,.08),97,k)
  old <- legacy_kernel(beta,se,engine)
  new <- apply_engine(beta,se,engine)
  for (field in names(old)) check(identical(old[[field]],new[[field]]))
}
cat('Exact-zero optimisation regression checks passed:', checks, '\n')
