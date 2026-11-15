# Patch for R SVN

[`stats-pacf-wei-trunk.patch`](stats-pacf-wei-trunk.patch) integrates the current
`pacfwei` 0.3.0 implementation into R SVN trunk **revision 90875**. It includes the C
and R implementations, mathematical documentation, author credit and email,
regression tests, and the Pinkham data and help with the Palda (1965) reference.

Apply it from the root of a full R source checkout at the matching revision:

```sh
git apply --check --whitespace=error-all /path/to/stats-pacf-wei-trunk.patch
git apply /path/to/stats-pacf-wei-trunk.patch
```

`git apply` also works in an SVN checkout. The pinned source revision is
<https://svn.r-project.org/R/!svn/rvr/90875/trunk/>. For a later revision, repeat
the application check and review any changes in the affected files.

After rebuilding R, the interface is:

```r
pacf(x, lag.max = 10, plot = FALSE)               # C implementation
pacf(x, lag.max = 10, plot = FALSE, method = "R")
```

The patch modifies `stats::pacf()` and adds Pinkham to `datasets`. It changes
the multivariate result from partial autoregression coefficients to Wei's
componentwise partial autocorrelations. `ar.yw(x)$partialacf` retains its
regression-coefficient interpretation. The univariate calculation retains
its existing algorithm.

Version 0.3.0 adds data-based SVD whitening, Cholesky solves, rank and residual
covariance safeguards, and safe binary prescaling. Residual matrices are mapped
back before normalization, preserving the estimator. Shared lag validation
accepts `Inf`, rejects malformed values consistently, and gives the default lag
a minimum of one. See [validation and compatibility notes](VALIDATION.md).

The patch passes strict application checks, standalone compilation, and
independent raw-data QR comparisons. A complete Windows build of R and its
base packages at revision 90875 succeeded; all 197 reviewer checks also passed
inside the rebuilt `stats`/`datasets` namespaces. Both help files pass Rd, HTML,
and plain-text conversion checks. The patch includes native registration and
Unix/Windows build entries. Other platforms and the full recommended-package
check suite remain outside this validation.

