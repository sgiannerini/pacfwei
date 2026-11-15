# Original algorithm: Copyright (c) 2007 Simone Giannerini; GPL >= 2.
# The original copyright and permission notice is retained in src/pacf_wei.c.
# The acf-class wrapper follows stats::pacf.default, Copyright (C) R Core Team,
# distributed under the GNU General Public License, version 2 or later.

.pacf_lag_max <- function(lag.max, n, m) {
    if (is.null(lag.max))
        lag.max <- max(1, floor(10 * log10(n / m)))
    if (!is.numeric(lag.max) || is.complex(lag.max) ||
        length(lag.max) != 1L || is.na(lag.max) || lag.max < 1)
        stop("'lag.max' must be a real numeric scalar of at least 1")
    as.integer(min(floor(lag.max), n - 1L))
}

.pacf_method_dots <- function(dots) {
    nm <- names(dots)
    if (is.null(nm)) return(invisible(NULL))
    partial <- !is.na(nm) & nzchar(nm) & nchar(nm) < nchar("method") &
        startsWith("method", nm)
    if (any(partial))
        stop(sprintf("argument '%s' must be named 'method' in full", nm[which(partial)[1L]]),
             call. = FALSE)
    invisible(NULL)
}

.pacf_prepare <- function(x, lag.max = NULL) {
    x <- as.matrix(x)
    if (!is.numeric(x) || is.complex(x) || any(!is.finite(x)))
        stop("'x' must contain finite real numeric values")
    n <- nrow(x)
    m <- ncol(x)
    if (n < 2L || m < 1L)
        stop("'x' must have at least two rows and one column")
    lag.max <- .pacf_lag_max(lag.max, n, m)
    snames <- colnames(x)
    magnitude <- apply(abs(x), 2L, max)
    if (any(magnitude == 0)) stop("zero-variance series in 'x'")
    # Binary scaling protects the range without rounding away large-offset
    # differences through an arbitrary divisor before centering.
    emin <- .Machine$double.min.exp - .Machine$double.digits + 1L
    emax <- .Machine$double.max.exp - 1L
    exponent <- pmax(emin, pmin(emax, floor(log2(magnitude))))
    x <- sweep(x, 2L, 2^exponent, "/")
    x <- sweep(x, 2L, colMeans(x), "-")
    magnitude <- apply(abs(x), 2L, max)
    if (any(!is.finite(magnitude) | magnitude <= 0))
        stop("zero-variance series in 'x'")
    x <- sweep(x, 2L, magnitude, "/")
    scales <- sqrt(colMeans(x * x))
    if (any(!is.finite(scales) | scales <= 0))
        stop("zero-variance series in 'x'")
    x <- sweep(x, 2L, scales, "/")
    map <- diag(m)
    if (lag.max > 1L) {
        if (m >= n) stop("numerically rank-deficient series in 'x'")
        decomposition <- svd(x, nu = m, nv = m)
        tolerance <- max(n, m) * .Machine$double.eps * decomposition$d[1L]
        if (any(!is.finite(decomposition$d)) ||
            any(decomposition$d <= tolerance))
            stop("numerically rank-deficient series in 'x'")
        # X = W %*% map. Covariances are estimated directly from W, rather
        # than whitening an already ill-conditioned covariance sequence.
        map <- sweep(t(decomposition$v), 1L, decomposition$d / sqrt(n), "*")
        x <- sqrt(n) * decomposition$u
    }
    list(cov = stats::acf(x, type = "covariance", lag.max = lag.max,
                         plot = FALSE, demean = FALSE)$acf,
         map = map, lag.max = lag.max, nser = m, snames = snames)
}

.pacf_R_kernel <- function(cov, map = NULL) {
    d <- dim(cov)
    if (!is.double(cov) || length(d) != 3L || d[1L] < 2L ||
        d[2L] < 1L || d[2L] != d[3L] || any(!is.finite(cov)))
        stop("'cov' must be a finite double array with dimensions (p+1,m,m)")
    p <- d[1L] - 1L
    m <- d[2L]
    if (is.null(map)) map <- diag(m)
    if (!is.double(map) || !identical(dim(map), c(m, m)) || any(!is.finite(map)))
        stop("'map' must be a finite double matrix with dimensions (m,m)")
    G <- lapply(seq_len(p + 1L), function(k) matrix(cov[k, , ], m, m))
    symmetry.tolerance <- 64 * m * .Machine$double.eps * max(1, norm(G[[1L]], "1"))
    if (max(abs(G[[1L]] - t(G[[1L]]))) > symmetry.tolerance)
        stop("lag-zero covariance must be symmetric")
    G[[1L]] <- (G[[1L]] + t(G[[1L]])) / 2
    vf <- vb <- G[[1L]]
    forward <- backward <- list()
    out <- array(0, c(p, m, m))
    coefficient.tolerance <- 64 * m * .Machine$double.eps
    covariance.check <- function(U, V, K, s) {
        H <- rbind(cbind(U, K), cbind(t(K), V))
        if (any(!is.finite(H)))
            stop(sprintf("nonfinite residual covariance at lag %d", s))
        tolerance <- 64 * (2 * m) * .Machine$double.eps * max(1, norm(H, "1"))
        if (min(eigen(H, symmetric = TRUE, only.values = TRUE)$values) < -tolerance)
            stop(sprintf("invalid joint residual covariance at lag %d", s))
    }
    factor <- function(A, s, side) {
        R <- tryCatch(chol(A), error = function(e)
            stop(sprintf("non-positive-definite %s residual covariance at lag %d",
                         side, s), call. = FALSE))
        condition <- rcond(A)
        if (!is.finite(condition) || condition < sqrt(.Machine$double.eps))
            stop(sprintf("numerically singular %s residual covariance at lag %d",
                         side, s), call. = FALSE)
        R
    }
    right.solve <- function(B, R) {
        t(backsolve(R, forwardsolve(t(R), t(B))))
    }
    for (s in seq_len(p)) {
        K <- G[[s + 1L]]
        if (s > 1L) for (k in seq_len(s - 1L))
            K <- K - forward[[k]] %*% G[[s - k + 1L]]
        covariance.check(vf, vb, K, s)
        U <- crossprod(map, vf %*% map)
        V <- crossprod(map, vb %*% map)
        H <- crossprod(map, K %*% map)
        u <- diag(U)
        v <- diag(V)
        if (any(!is.finite(c(u, v))) || any(c(u, v) <= 0))
            stop(sprintf("nonpositive residual variance at lag %d", s))
        value <- H / outer(sqrt(u), sqrt(v))
        if (any(!is.finite(value)))
            stop(sprintf("nonfinite partial correlation at lag %d", s))
        if (any(abs(value) > 1 + coefficient.tolerance))
            stop(sprintf("partial correlation exceeds its numerical bound at lag %d", s))
        out[s, , ] <- pmax(-1, pmin(1, value))
        if (s == p) break
        Rf <- factor(vf, s, "forward")
        Rb <- factor(vb, s, "backward")
        Fs <- right.solve(K, Rb)
        Bs <- right.solve(t(K), Rf)
        old.forward <- forward
        old.backward <- backward
        if (s > 1L) for (k in seq_len(s - 1L)) {
            forward[[k]] <- old.forward[[k]] - Fs %*% old.backward[[s-k]]
            backward[[k]] <- old.backward[[k]] - Bs %*% old.forward[[s-k]]
        }
        forward[[s]] <- Fs
        backward[[s]] <- Bs
        vf <- vf - Fs %*% t(K)
        vb <- vb - Bs %*% K
        vf <- (vf + t(vf)) / 2
        vb <- (vb + t(vb)) / 2
    }
    out
}

.pacf_mts <- function(x, lag.max, method) {
    z <- .pacf_prepare(x, lag.max)
    out <- if (method == "C") .Call(C_pacf_wei, z$cov, z$map) else
        .pacf_R_kernel(z$cov, z$map)
    dimnames(out) <- list(as.character(seq_len(z$lag.max)), z$snames, z$snames)
    out
}

#' Componentwise partial autocorrelations for vector time series
#'
#' Computes Wei's partial autocorrelations using a stabilized forward/backward
#' recursion. Select the native C implementation (the default) or its R reference.
#' Results use the lag convention and plotting interface of `stats::pacf()`.
#'
#' @param x A real numeric vector, matrix, data frame or time series. After
#'   `na.action`, at least two observations are required. Multivariate data must
#'   be finite, complete and have positive marginal variances.
#' @param lag.max Maximum lag: `NULL`, or a real numeric scalar of at least one.
#'   The default is `max(1, floor(10 * log10(n / m)))`, where \eqn{n}{n} is the
#'   number of observations and \eqn{m}{m} the number of components. Positive
#'   fractional values are floored and capped at \eqn{n-1}{n - 1}; `Inf` requests
#'   that cap. The same validation applies to univariate and multivariate input.
#'   Capping does not guarantee that every requested lag is numerically estimable.
#' @param plot Whether to plot using the standard `acf` method.
#' @param na.action Missing-value handler, applied before calculation.
#'   Passing remaining missing values is not supported for multivariate data.
#' @param ... Further arguments to the `acf` plotting method. Abbreviations of
#'   `method`, such as `meth`, are rejected without evaluating their expressions.
#' @param method Multivariate recursion: `"C"` (default) or `"R"`. This argument
#'   must be named in full because it follows `...`. Univariate input uses
#'   `stats::pacf()` for either choice, after the shared argument validation.
#'
#' @details
#' \strong{Definition.} Let \eqn{Z_t}{Z[t]} be a centered stationary vector process.
#' At lag \eqn{s}{s}, project the observations at times \eqn{t+s}{t+s} and
#' \eqn{t}{t} on every intervening vector and form residuals
#' \deqn{u_{s-1,t+s}=Z_{t+s}-\sum_{k=1}^{s-1}\alpha_{s-1,k}Z_{t+s-k}.}{
#' u = Z[t+s] - sum(k = 1,...,s-1) alpha[s-1,k] Z[t+s-k].}
#' \deqn{v_{s-1,t}=Z_t-\sum_{k=1}^{s-1}\beta_{s-1,k}Z_{t+k}.}{
#' v = Z[t] - sum(k = 1,...,s-1) beta[s-1,k] Z[t+k].}
#' These are linear projections, without a Gaussian distribution assumption.
#' Write \eqn{U_s=\operatorname{Var}(u)}{U[s] = Var(u)},
#' \eqn{V_s=\operatorname{Var}(v)}{V[s] = Var(v)} and
#' \eqn{K_s=\operatorname{Cov}(u,v)}{K[s] = Cov(u, v)}. Let
#' \eqn{D_U(s)}{D_U[s]} and \eqn{D_V(s)}{D_V[s]} be diagonal matrices of
#' marginal residual standard deviations. The returned lag matrix is
#' \deqn{P_R(s)=D_U(s)^{-1}K_sD_V(s)^{-1}.}{
#' P_R[s] = inverse(D_U[s]) K[s] inverse(D_V[s]).}
#' Entry \eqn{(i,j)}{(i,j)} correlates the residual for component \eqn{i}{i}
#' at time \eqn{t+s}{t+s} with that for component \eqn{j}{j} at time \eqn{t}{t}.
#' There is no additional conditioning on the other components at those times.
#' The lag matrices are generally nonsymmetric.
#'
#' \strong{Relation to Wei.} Wei (2006, Section 16.5.3, pp. 408--414) uses
#' \eqn{\Gamma(s)=\operatorname{Cov}(Z_t,Z_{t+s})}{Gamma[s] = Cov(Z[t], Z[t+s])}
#' and the transpose of the returned matrix in (16.5.43):
#' \eqn{P_R(s)=P(s)'}{P_R[s] = transpose(P[s])}.
#' This follows R's `stats::acf()` convention: `result$acf[s + 1, i, j]` uses
#' component \eqn{i}{i} at time \eqn{t+s}{t+s} and component \eqn{j}{j} at time
#' \eqn{t}{t}. The partial autoregression matrix in (16.5.44) is instead
#' \eqn{K_s V_s^{-1}}{K[s] inverse(V[s])}; its entries need not lie in
#' \eqn{[-1,1]}{[-1,1]}. The definitions coincide for a univariate process.
#'
#' \strong{Recursion.} With
#' \eqn{C_h=\operatorname{Cov}(Z_{t+h},Z_t)}{C[h] = Cov(Z[t+h], Z[t])},
#' initialize \eqn{U_1=V_1=C_0}{U[1] = V[1] = C[0]}. Using coefficient
#' matrices from the previous order, compute
#' \deqn{K_s=C_s-\sum_{k=1}^{s-1}F_{s-1,k}C_{s-k}.}{
#' K[s] = C[s] - sum(k = 1,...,s-1) F[s-1,k] C[s-k].}
#' \deqn{F_{s,s}=K_s V_s^{-1},\quad B_{s,s}=K_s' U_s^{-1}.}{
#' F[s,s] = K[s] inverse(V[s]); B[s,s] = transpose(K[s]) inverse(U[s]).}
#' For \eqn{1\leq k<s}{1 <= k < s}, update both sets from the old order:
#' \deqn{F_{s,k}=F_{s-1,k}-F_{s,s}B_{s-1,s-k}.}{
#' F[s,k] = F[s-1,k] - F[s,s] B[s-1,s-k].}
#' \deqn{B_{s,k}=B_{s-1,k}-B_{s,s}F_{s-1,s-k}.}{
#' B[s,k] = B[s-1,k] - B[s,s] F[s-1,s-k].}
#' The residual covariances for the next order are
#' \deqn{U_{s+1}=U_s-F_{s,s}K_s',\quad V_{s+1}=V_s-B_{s,s}K_s.}{
#' U[s+1] = U[s] - F[s,s] transpose(K[s]); V[s+1] = V[s] - B[s,s] K[s].}
#' These are (16.5.50)--(16.5.51) in R's covariance convention.
#' The displayed inverses specify the mathematics; both implementations solve
#' systems using Cholesky factors without explicitly forming matrix inverses.
#'
#' \strong{Numerical stabilization.} Power-of-two prescaling protects the
#' floating-point range before centering. Arbitrary marginal normalization
#' follows centering, avoiding loss of large-offset differences. Covariance
#' estimates use `stats::acf()` with denominator \eqn{n}{n}.
#' For more than one lag, a full-rank SVD of the centered, marginally normalized
#' observation matrix provides \eqn{X=W A}{X = W A}, with
#' \eqn{W'W/n=I}{transpose(W) W / n = I}. The recursion uses covariances
#' estimated directly from \eqn{W}{W}. Before normalization, its residual
#' matrices are transformed back:
#' \deqn{U_X=A'U_W A,\quad V_X=A'V_W A,\quad K_X=A'K_W A.}{
#' U_X = transpose(A) U_W A; V_X = transpose(A) V_W A; K_X = transpose(A) K_W A.}
#' An invertible transformation preserves the span of all intervening vector
#' observations, so this leaves the estimator unchanged in exact arithmetic.
#' Transforming already normalized correlations would not be equivalent.
#'
#' No singular directions are silently discarded and no regularization is
#' added. SVD singular values must exceed
#' \eqn{\max(n,m)\epsilon d_{\max}}{max(n,m) * eps * largest_singular_value},
#' where \eqn{\epsilon}{eps} is `.Machine$double.eps`. Residual systems that
#' require a solve must be positive definite with estimated reciprocal condition
#' number at least \eqn{\sqrt{\epsilon}}{sqrt(eps)} in whitened coordinates.
#' Invalid joint residual covariances are rejected, including at the final lag;
#' eigenvalues down to
#' \eqn{-64(2m)\epsilon\max(1,\|H\|_1)}{-64 * (2*m) * eps * max(1, norm(H, "1"))}
#' are treated as roundoff for the joint residual block \eqn{H}{H}.
#' Exact positive-semidefinite boundaries are allowed. At lag one no SVD or
#' regression solve is needed, so redundant columns with positive variances
#' are permitted. A coefficient exceeding one in modulus by at most
#' \eqn{64m\epsilon}{64 * m * eps} is clamped to the boundary; a larger
#' violation raises an error. Short or highly predictable series can still
#' produce singular or numerically unreliable higher orders.
#'
#' Dense arithmetic requires \eqn{O(nm^2+m^3p^2)}{O(n*m^2 + m^3*p^2)} operations
#' and \eqn{O(nm+m^2p)}{O(n*m + m^2*p)} storage, where \eqn{p}{p} is the
#' maximum lag. C and R can differ by roundoff and borderline failure decisions.
#'
#' This multivariate definition differs from the partial autoregression
#' coefficients returned by `stats::pacf()` and `stats::ar.yw(x)$partialacf`
#' in the R version for which this package was developed.
#' Loading this package changes neither function.
#'
#' @return An object of class `"acf"`, with components `acf`, `lag`, `type`,
#'   `n.used`, `series` and `snames`. For a result named `result`, `result$acf`
#'   is a \eqn{p\times m\times m}{p x m x m} coefficient array, ordered by lag,
#'   component at time \eqn{t+s}{t+s}, and component at time \eqn{t}{t}.
#'   Lags are in time units; plotting uses the standard signed cross-lag
#'   convention. Returned invisibly if `plot = TRUE`.
#' @author Simone Giannerini \email{simone.giannerini@uniud.it}
#' @references
#' Wei, W. W. S. (2006). \emph{Time Series Analysis: Univariate and Multivariate
#' Methods}, second edition. Pearson Prentice Hall. Section 16.5.3,
#' equations (16.5.39)--(16.5.44) and (16.5.50)--(16.5.51).
#' @seealso \code{\link[stats]{acf}}, \code{\link[stats]{pacf}},
#'   \code{\link[stats]{ar.yw}}, \code{\link{pinkham}}
#' @examples
#' set.seed(10)
#' x <- cbind(x = rnorm(1000, sd = 10000), y = rnorm(1000))
#' pc <- pacf_wei(x, lag.max = 10, plot = FALSE)
#' pr <- pacf_wei(x, lag.max = 10, plot = FALSE, method = "R")
#' stopifnot(isTRUE(all.equal(pc$acf, pr$acf, tolerance = 1e-10)))
#' aperm(pc$acf, c(1, 3, 2))  # Wei's theoretical orientation
#' @export
#' @md
#' @useDynLib pacfwei, .registration = TRUE, .fixes = "C_"

pacf_wei <- function(x, lag.max = NULL, plot = TRUE,
                     na.action = stats::na.fail, ..., method = c("C", "R")) {
    .pacf_method_dots(match.call(expand.dots = FALSE)[["..."]])
    method <- match.arg(method)
    series <- deparse1(substitute(x))
    x <- na.action(stats::as.ts(x))
    if (!is.numeric(x) || is.complex(x)) stop("'x' must be real numeric")
    if (NROW(x) < 2L) stop("'x' must have at least two observations")
    x <- drop(x)
    lag.max <- .pacf_lag_max(lag.max, NROW(x), NCOL(x))
    if (!is.matrix(x)) {
        out <- stats::pacf(x, lag.max = lag.max, plot = FALSE,
                           na.action = stats::na.pass)
        out$series <- series
    } else {
        values <- .pacf_mts(x, lag.max, method)
        p <- dim(values)[1L]
        m <- ncol(x)
        signs <- matrix(1, m, m)
        signs[lower.tri(signs)] <- -1
        out <- structure(list(acf = values, type = "partial", n.used = nrow(x),
                              lag = outer(seq_len(p), signs / stats::frequency(x)),
                              series = series, snames = colnames(x)), class = "acf")
    }
    if (plot) {
        graphics::plot(out, ...)
        invisible(out)
    } else out
}
