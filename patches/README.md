# Patch for R SVN

[`stats-pacf-wei-trunk.patch`](stats-pacf-wei-trunk.patch) integrates the current
`pacfwei` implementation into R SVN trunk **revision 90875**. It includes the C
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

The patch passed strict application checks, standalone C compilation, and
isolated numerical tests, including C/R agreement, independent projection
calculations, and exact agreement with the package. The data retain the
package's values and data-frame attributes. Both help files passed Rd and
HTML conversion checks.

**A full build and check of R itself has not been performed.** Rebuild R,
including `stats` and `datasets`, and run its checks before proposing upstream
integration. The patch includes native registration and Unix/Windows build
entries. It is a proposed change and has not been accepted into R.

The alternative patch for the supplied local source archive is not retained.
