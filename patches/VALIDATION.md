# pacfwei 0.3.0: numerical revision

This revision addresses the numerical and interface feedback on the initial
R SVN patch. C remains the default; the public package function remains
`pacf_wei()`.

## Changes and numerical policy

The recursion now operates on covariances estimated directly from SVD-whitened
observations. Residual covariance and cross-covariance matrices are mapped back
before computing componentwise correlations. This preserves the estimator in
exact arithmetic while avoiding the ill-conditioned covariance representation.

Both implementations use Cholesky solves, reject numerical rank deficiency
without dropping directions, and validate the full joint residual covariance
at every lag. The final lag may have singular marginal residual matrices if
the requested componentwise correlations are still defined. Lag one requires
no whitening or regression solve and continues to support redundant columns.

For reproducibility, the numerical policies are explicit:

- SVD rank threshold: `max(n, m) * eps * largest_singular_value`.
- Minimum reciprocal condition estimate for systems solved in whitened
  coordinates: `sqrt(eps)`.
- Joint covariance eigenvalue tolerance:
  `64 * (2*m) * eps * max(1, norm(H, "1"))`.
- Correlation overshoot tolerance: `64 * m * eps`. Only such tiny overshoots
  are clamped to -1 or 1; larger violations raise an error.

Power-of-two prescaling precedes centering; arbitrary magnitude and RMS
normalization follow it. Lag validation is shared across uni/multivariate
inputs, accepts `Inf`, rejects malformed values consistently, and floors the
default at one. These stricter argument rules are explicit compatibility
changes. Abbreviated method names are diagnosed without evaluating plotting
arguments. Plain-text help has ASCII mathematical alternatives.

## Accuracy checks on Windows / R 4.6.0

- `R CMD check --no-manual --no-vignettes`: **Status: OK**.
- Ordinary mathematical, invariance, metadata, native-input, and Pinkham suite:
  **66 checks for each implementation**.
- Reviewer regression suite: **197 checks**, including 25 nearly collinear
  samples, invalid joint covariances, exact semidefinite boundaries, large
  offsets, huge/subnormal scales, and shared argument handling.
- Independent raw-data QR oracle: largest difference **1.14e-10** across the
  25 near-collinear samples; largest C/R difference **5.55e-17**.
- Additional 24-setting grid: sample lengths 50/200, dimensions 2/4, lags 5/10,
  and collinearity noise 1e-3/1e-6/1e-8. Largest QR difference **2.00e-8**,
  largest C/R difference **1.11e-16**. Raw-data QR is itself increasingly
  sensitive as the collinearity noise shrinks.
- The reviewer examples now return max absolute coefficients **0.27691528222**
  for seed 4 and **0.165914811611** for seed 83.

The QR oracle uses zero-padded lagged raw observations, whose Gram matrix equals
the biased sample block-Toeplitz covariance used by `acf()`. It avoids both
the recurrence and an ill-conditioned covariance inversion. Bounds and C/R
agreement alone would not establish accuracy.

The SVN patch also passes strict `git apply --check --whitespace=error-all`
and isolated integration checks against revision 90875. Exact Pinkham values
and attributes are preserved. Package help and both patched help files have
been checked as Rd, HTML, and plain text.

## Rebuilt R SVN revision 90875

The complete Windows core-R and base-package build (`make -j4 all`) completed
successfully with Rtools45. Source commit
`9e1decef2ba9b4335828def1154c017c5b01b1a5` in the R source mirror corresponds
to SVN trunk revision 90875; the patch applied strictly to that source tree.
The installed Tcl/Tk 8.6 development bundle supplied the optional Windows
Tcl/Tk build dependency.

The 197 reviewer checks also passed inside the rebuilt R 4.7.0 runtime using
its own `stats` and `datasets` namespaces. The largest raw-QR difference was
1.60e-10 and the largest C/R difference 1.67e-16. Additional checks verified
the two-argument native registration, generic `acf(type = "partial")`
dispatch, exact Pinkham data, and terminal mathematical help. This validates
integration beyond loading the C kernel in an isolated R 4.6 environment.

The base examples were exercised, with expected Pinkham output baselines
updated in the patch. Remaining example comparison notes concern unrelated
platform text wrapping, encoding, and unbuilt vignette metadata. Network tests,
recommended-package-dependent tests, `check-all`, and other operating systems
are not covered. The initial broad check stopped at missing `Matrix`; the
base-only continuation uses explicit finite test targets rather than an empty
make-target list. Logs are retained with the release artifacts.

The selected specific, encoding, and general regression tests passed
individually. Running the general regressions concurrently caused a shared
`myLib` fixture conflict in `reg-packages.R`; rerunning that script alone
completed successfully (exit status 0). Thus the parallel make invocation is
not reported as a clean aggregate check, and this validation is not a claim
that R's complete `check-all` suite passed.

## Local timings

Median elapsed milliseconds per complete call, seven randomized paired blocks
with warm-up and batching targeting approximately 0.3 seconds. Data are
independent Gaussian series; no compiler or other validation jobs ran during
this measurement. These are within-session comparisons, not estimates to
compare directly with older timing sessions.

| n | m | lags | Previous C (0.2.0) | Revised C | Revised R |
|---:|---:|---:|---:|---:|---:|
| 1000 | 2 | 10 | 0.99 | 1.35 | 3.16 |
| 2000 | 4 | 30 | 3.70 | 4.67 | 9.12 |
| 5000 | 10 | 50 | 48.57 | 55.71 | 65.00 |
| 5000 | 50 | 10 | 245.00 | 275.00 | 285.00 |

The SVD and covariance safeguards add work; C remains faster than the R
reference in these cases. Raw block timings and the additional accuracy grid
are cached as CSV files.
