# Numerical cases identified in the review, checked against raw-data projection.
# Usage: Rscript tests/reviewer-regressions.R [installed-library]
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && nzchar(args[1L])) .libPaths(c(args[1L], .libPaths()))
stopifnot(requireNamespace("pacfwei", quietly = TRUE))
pacf_wei <- pacfwei::pacf_wei
nchecks <- 0L
check <- function(ok, label) {
    if (!isTRUE(ok)) stop("FAILED: ", label, call. = FALSE)
    nchecks <<- nchecks + 1L
}
near <- function(a, b, label, tolerance = 2e-12) {
    err <- if (length(a) == length(b)) max(abs(as.numeric(a) - as.numeric(b))) else Inf
    check(identical(dim(a), dim(b)) && is.finite(err) && err <= tolerance,
          sprintf("%s (max error %.3g)", label, err))
    invisible(err)
}
throws <- function(expr, label, pattern = NULL) {
    value <- tryCatch({force(expr); NULL}, error = identity)
    check(inherits(value, "error") &&
          (is.null(pattern) || grepl(pattern, conditionMessage(value))), label)
}

# The Gram matrix of these zero-padded shifts is n times the biased sample
# block-Toeplitz covariance. Householder QR projects the two observations on
# every intervening vector without constructing or inverting that covariance.
# This oracle does not whiten or use the forward/backward recursion.
qr_reference <- function(x, p) {
    x <- as.matrix(x)
    x <- sweep(x, 2L, colMeans(x), "-")
    n <- nrow(x); m <- ncol(x)
    out <- array(0, c(p, m, m))
    for (s in seq_len(p)) {
        padded <- matrix(0, n + s, (s + 1L) * m)
        for (k in 0:s)
            padded[k + seq_len(n), k * m + seq_len(m)] <- x
        later <- padded[, seq_len(m), drop = FALSE]
        earlier <- padded[, s * m + seq_len(m), drop = FALSE]
        if (s > 1L) {
            middle <- padded[, seq.int(m + 1L, s * m), drop = FALSE]
            Q <- qr.Q(qr(middle, LAPACK = TRUE))
            later <- later - Q %*% crossprod(Q, later)
            earlier <- earlier - Q %*% crossprod(Q, earlier)
        }
        out[s, , ] <- crossprod(later, earlier) /
            outer(sqrt(colSums(later * later)), sqrt(colSums(earlier * earlier)))
    }
    out
}

worst.oracle <- worst.backend <- 0
# Includes the reported cases and all five Cholesky-only survivors identified
# during exploratory testing. A bound check alone would miss their old errors.
seeds <- unique(c(4L, 83L, 14L, 134L, 237L, 387L, 954L, seq_len(20L)))
for (seed in seeds) {
    set.seed(seed)
    v <- rnorm(50)
    x <- cbind(v, v + 1e-6 * rnorm(50))
    expected <- qr_reference(x, 10L)
    cval <- pacf_wei(x, 10L, plot = FALSE, method = "C")$acf
    rval <- pacf_wei(x, 10L, plot = FALSE, method = "R")$acf
    worst.oracle <- max(worst.oracle,
        near(cval, expected, paste("nearly collinear C vs QR, seed", seed), 1e-8),
        near(rval, expected, paste("nearly collinear R vs QR, seed", seed), 1e-8))
    worst.backend <- max(worst.backend,
        near(cval, rval, paste("nearly collinear C/R agreement, seed", seed)))
    check(max(abs(cval)) <= 1 && max(abs(rval)) <= 1,
          paste("nearly collinear coefficients are bounded, seed", seed))
}

z <- cbind(z1 = c(-4, -2, -1, 0, 1, 2, 4, -3, 3, 0),
           z2 = c(1, -1, 2, -2, 3, -3, 4, -4, 5, -5))
offset <- cbind(z1 = 5.7e15 + z[, 1L], z2 = z[, 2L])
stopifnot(identical(sweep(offset, 2L, colMeans(offset), "-"), z))
set.seed(9)
u <- rnorm(50); u <- u - mean(u)
perfect <- cbind(c(u, 0), c(0, u))
smallest <- .Machine$double.xmin * .Machine$double.eps
subnormal <- cbind(1:10, -(1:10)) * smallest
huge <- cbind(c(.Machine$double.xmax, rep(-.Machine$double.xmax, 9L)), z[, 2L])
huge.scaled <- huge
huge.scaled[, 1L] <- huge.scaled[, 1L] / 2^(.Machine$double.max.exp - 1L)
set.seed(209)
regular <- matrix(rnorm(80), 40L, 2L)
wide <- matrix(rnorm(80), 10L, 8L)
for (method in c("C", "R")) {
    calc <- function(x, lag.max = NULL, ...) pacf_wei(x, lag.max,
        plot = FALSE, method = method, ...)$acf
    near(calc(offset, 1L)[1L, , ], stats::acf(offset, 1L, plot = FALSE)$acf[2L, , ],
         paste(method, "large-offset lag one matches ordinary ACF"), 1e-14)
    near(calc(offset, 1L), calc(z, 1L), paste(method, "exact-offset invariance"), 1e-14)
    check(max(abs(calc(perfect, 1L))) <= 1, paste(method, "perfect lag-one bound"))
    near(calc(subnormal, 1L), calc(cbind(1:10, -(1:10)), 1L),
         paste(method, "subnormal scaling invariance"), 1e-14)
    near(calc(huge, 1L), calc(huge.scaled, 1L),
         paste(method, "finite extreme input avoids centering overflow"), 1e-14)
    check(all(is.finite(calc(huge, 1L))), paste(method, "finite huge-scale result"))
    check(identical(dim(calc(wide)), c(1L, 8L, 8L)),
          paste(method, "short-wide default is lag one"))
    for (input in list(regular, regular[, 1L], regular[, 1L, drop = FALSE])) {
        near(calc(input, Inf), calc(input, nrow(regular) - 1L),
             paste(method, "Inf caps identically in uni/multivariate input"))
        near(calc(input, 3.9), calc(input, 3L), paste(method, "fractional lag is floored"))
        for (bad in list(NA_real_, NaN, -Inf, 0, -1, "5", c(2, 5), TRUE, numeric()))
            throws(calc(input, bad), paste(method, "malformed lag rejected consistently"), "lag.max")
    }
    throws(pacf_wei("a", plot = FALSE, method = method),
           paste(method, "numeric-type validation precedes length check"), "numeric")
    throws(pacf_wei(matrix(1:3, 1L, 3L), plot = FALSE, method = method),
           paste(method, "one-observation matrix rejected before drop"), "observations")
}
for (plot in c(FALSE, TRUE)) {
    throws(pacf_wei(regular, 2L, plot = plot, meth = "invalid"),
           "abbreviated method rejected even when plotting is disabled", "method")
    throws(pacf_wei(regular, 2L, plot = plot, meth = stop("must not evaluate")),
           "method-prefix validation inspects unevaluated dots", "method")
}

native.symbol <- get("C_pacf_wei", asNamespace("pacfwei"))
native <- function(cov, map = diag(dim(cov)[2L])) .Call(native.symbol, cov, map)
joint <- array(0, c(2L, 2L, 2L))
joint[1L, , ] <- diag(2L)
# U and V are individually SPD and every entry of K is below one, yet the
# full joint residual covariance is indefinite (largest singular value 1.5).
joint[2L, , ] <- matrix(0.75, 2L, 2L)
throws(native(joint), "invalid final-lag joint covariance fails")
joint[2L, , ] <- diag(2L)
near(native(joint)[1L, , ], diag(2L), "semidefinite exact-boundary joint block accepted")
duplicate <- array(0, c(2L, 2L, 2L))
duplicate[1L, , ] <- matrix(1, 2L, 2L)
duplicate[2L, , ] <- matrix(0.2, 2L, 2L)
check(all(is.finite(native(duplicate))), "singular marginal covariance permitted at lag one")
boundary <- array(0, c(3L, 2L, 2L))
boundary[1L, , ] <- diag(2L)
boundary[2L, , ] <- matrix(c(0, 0, 1, 0), 2L, 2L)
rotation <- matrix(c(1, -1, 1, 1), 2L, 2L) / sqrt(2)
# At the final lag, each whitened marginal residual matrix is singular,
# while the original components retain positive variances after mapping.
near(native(boundary, rotation)[2L, , ], matrix(0, 2L, 2L),
     "semidefinite final marginal covariances accepted when correlations exist")
near(getFromNamespace(".pacf_R_kernel", "pacfwei")(boundary, rotation),
     native(boundary, rotation), "R/C agree at a singular final residual boundary")
throws(native(joint, matrix(1, 2L, 3L)), "native rejects wrong map dimensions")
throws(native(joint, matrix(1L, 2L, 2L)), "native rejects integer map")
bad.map <- diag(2L); bad.map[1L] <- NA_real_
throws(native(joint, bad.map), "native rejects nonfinite map")
throws(native(joint, numeric(4L)), "native rejects undimensioned map")

cat(sprintf("All %d reviewer regression checks passed.\n", nchecks))
cat(sprintf("Nearly collinear grid: %d seeds; max QR error %.3g; max C/R difference %.3g.\n",
            length(seeds), worst.oracle, worst.backend))
