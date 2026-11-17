# SubtractionSEM

Fast GWAS-by-subtraction for square, invertible first-order factor models.
Uses the fixed-loading and sampling-covariance conventions of analytic GenomicSEM.

## Install

Requires R >= 4.1 and the usual R build tools. Install the validated dependencies:

```r
install.packages("remotes")
remotes::install_github("cran/lavaan@8c44109a2856a08299e1ca71b409c3c6f2d2195a", upgrade = "never")
remotes::install_github("GenomicSEM/GenomicSEM@6b65ca5db39fdade08b0d811477be1cdd57b5039", upgrade = "never")
remotes::install_github("gpmerola/SubtractionSEM", upgrade = "never")
```

Keep lavaan at 0.6-21 for compatibility with the reference fitter.

## Use

```r
library(SubtractionSEM)
result <- subtractionGWAS(covstruc, SNPs, model, batch_size = 100000)
# Optional: reuse fitted loadings with usermod = fitted_model
?subtractionGWAS
```

`covstruc` is a joint GenomicSEM LDSC result containing S, V and I.
`SNPs` is harmonised `GenomicSEM::sumstats()` output with `beta.<trait>` and
`se.<trait>` columns matching S. `model` uses lavaan syntax with as many
first-order factors as traits and one free `~ SNP` path per factor.

Results preserve SNP order and metadata, with `beta_`, `SE_`, `Z_beta_` and
`p_val_` columns for each factor. Invalid input rows retain their positions
with missing estimates and a `subtraction_status` explanation.

Loadings are treated as fixed: their estimation uncertainty is not propagated.
Alleles must already be harmonised. No omnibus Q statistic is returned.
Use GenomicSEM for models outside the supported class. See the R help page for
input requirements and numerical safeguards.

## Check

With dependencies installed, run from the repository root:

```sh
R CMD build --no-build-vignettes .
R CMD check --no-manual SubtractionSEM_0.2.0.tar.gz
```

The synthetic tests compare results with the pinned GenomicSEM reference.
No empirical data are included. GitHub Actions runs these checks on Linux,
macOS and Windows.

## License

Copyright (C) 2026 Giuseppe Pierpaolo Merola. GPL version 3; see LICENSE.
The code, documentation and tests may be redistributed and modified under
that license, without warranty, including merchantability or fitness for a
particular purpose. GenomicSEM retains its upstream GPL licensing and attribution.
SubtractionSEM is maintained separately and does not imply upstream endorsement.
