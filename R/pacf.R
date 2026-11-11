# Original algorithm: Copyright (c) 2007 Simone Giannerini; GPL >= 2.
# The original copyright and permission notice is retained in src/pacf_wei.c.
# The acf-class wrapper follows stats::pacf.default, Copyright (C) R Core Team,
# distributed under the GNU General Public License, version 2 or later.

.pacf_prepare <- function(x, lag.max = NULL) {
    x <- as.matrix(x)
    if (!is.numeric(x) || is.complex(x) || any(!is.finite(x)))
        stop("'x' must contain finite real numeric values")
    n <- nrow(x)
    m <- ncol(x)
    if (n < 2L || m < 1L)
        stop("'x' must have at least two rows and one column")
    if (is.null(lag.max)) lag.max <- floor(10 * log10(n / m))
    if (length(lag.max) != 1L || !is.numeric(lag.max) ||
        is.complex(lag.max) || !is.finite(lag.max) || lag.max < 1)
        stop("'lag.max' must be a finite number of at least 1")
    lag.max <- as.integer(min(floor(lag.max), n - 1L))
    magnitude <- apply(abs(x), 2L, max)
    if (any(magnitude == 0)) stop("zero-variance series in 'x'")
    x <- sweep(x, 2L, magnitude, "/")
    x <- sweep(x, 2L, colMeans(x), "-")
    scales <- sqrt(colMeans(x * x))
    if (any(!is.finite(scales) | scales <= 0))
        stop("zero-variance series in 'x'")
    x <- sweep(x, 2L, scales, "/")
    list(cov = stats::acf(x, type = "covariance", lag.max = lag.max,
                         plot = FALSE, demean = FALSE)$acf,
         lag.max = lag.max, nser = m, snames = colnames(x))
}

.pacf_R_kernel <- function(cov) {
    d <- dim(cov)
    if (!is.double(cov) || length(d) != 3L || d[1L] < 2L ||
        d[2L] < 1L || d[2L] != d[3L] || any(!is.finite(cov)))
        stop("'cov' must be a finite double array with dimensions (p+1,m,m)")
    p <- d[1L] - 1L
    m <- d[2L]
    G <- lapply(seq_len(p + 1L), function(k) matrix(cov[k, , ], m, m))
    vf <- vb <- G[[1L]]
    forward <- backward <- list()
    out <- array(0, c(p, m, m))
    right.solve <- function(B, A, s) {
        tryCatch(t(solve(t(A), t(B))), error = function(e)
            stop(sprintf("singular residual covariance at lag %d: %s",
                         s, conditionMessage(e)), call. = FALSE))
    }
    for (s in seq_len(p)) {
        K <- G[[s + 1L]]
        if (s > 1L) for (k in seq_len(s - 1L))
            K <- K - forward[[k]] %*% G[[s - k + 1L]]
        u <- diag(vf)
        v <- diag(vb)
        if (any(!is.finite(c(u, v))) || any(c(u, v) <= 0))
            stop(sprintf("nonpositive residual variance at lag %d", s))
        value <- K / outer(sqrt(u), sqrt(v))
        if (any(!is.finite(value)))
            stop(sprintf("nonfinite partial correlation at lag %d", s))
        out[s, , ] <- value
        if (s == p) break
        Fs <- right.solve(K, vb, s)
        Bs <- right.solve(t(K), vf, s)
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
    out <- if (method == "C") .Call(C_pacf_wei, z$cov) else .pacf_R_kernel(z$cov)
    dimnames(out) <- list(as.character(seq_len(z$lag.max)), z$snames, z$snames)
    out
}

#' Componentwise partial autocorrelations for vector time series
#'
#' Computes Wei's partial autocorrelations using either native C or an R
#' implementation of the same recursion. The C implementation is the default.
#' Results use the lag convention and plotting interface of `stats::pacf()`.
#'
#' @param x A real numeric vector, matrix, data frame or time series. After
#'   `na.action`, at least two observations are required. Multivariate data must
#'   be finite, complete and have positive marginal variances.
#' @param lag.max Maximum lag. For multivariate input, defaults to
#'   `floor(10 * log10(n / m))`, is floored and capped at \eqn{n-1}, and must be
#'   finite and at least one. Here \eqn{n} is the observation count and \eqn{m}
#'   is the number of components. Univariate lag handling follows `stats::pacf()`.
#' @param plot Whether to plot using the standard `acf` method.
#' @param na.action Missing-value handler, applied before calculation.
#'   Passing remaining missing values is not supported for multivariate data.
#' @param ... Further arguments to the `acf` plotting method.
#' @param method Multivariate recursion: `"C"` (default) or `"R"`. This argument
#'   must be named because it follows `...`. Univariate input uses the existing
#'   `stats::pacf()` implementation for either choice.
#'
#' @details
#' \strong{Definition.} Let \eqn{Z_t} be a centered stationary vector process.
#' At lag \eqn{s}, project the vector observations \eqn{Z_{t+s}} and \eqn{Z_t}
#' on all the intervening vectors \eqn{Z_{t+1},\ldots,Z_{t+s-1}} and form residuals
#' \deqn{u_{s-1,t+s}=Z_{t+s}-\sum_{k=1}^{s-1}\alpha_{s-1,k}Z_{t+s-k},}
#' \deqn{v_{s-1,t}=Z_t-\sum_{k=1}^{s-1}\beta_{s-1,k}Z_{t+k}.}
#' These are linear projections, not an assumption of Gaussian conditional
#' distributions. Write \eqn{U_s=\operatorname{Var}(u_{s-1,t+s})},
#' \eqn{V_s=\operatorname{Var}(v_{s-1,t})} and
#' \eqn{K_s=\operatorname{Cov}(u_{s-1,t+s},v_{s-1,t})}.
#' Let \eqn{D_U(s)} and \eqn{D_V(s)} be diagonal matrices of marginal residual
#' standard deviations. In R's lag convention, the returned matrix is
#' \deqn{P_R(s)=D_U(s)^{-1}K_sD_V(s)^{-1}.}
#' Consequently entry \eqn{(i,j)} is the correlation between the residuals for
#' component \eqn{i} at time \eqn{t+s} and component \eqn{j} at time \eqn{t},
#' after projecting out the intervening vectors. The individual lag matrix is
#' generally nonsymmetric. There is no additional conditioning on the other
#' components at time \eqn{t} or \eqn{t+s}.
#'
#' \strong{Relation to Wei.} Wei (2006, Section 16.5.3, pp. 408--414) uses
#' \eqn{\Gamma(s)=\operatorname{Cov}(Z_t,Z_{t+s})} and defines
#' \eqn{P(s)=D_V(s)^{-1}V_{vu}(s)D_U(s)^{-1}} in (16.5.43),
#' with \eqn{V_{vu}(s)=K_s'}.
#' Thus \eqn{P_R(s)=P(s)'}. This transpose is required for consistency with
#' `stats::acf()`, whose positive-lag entry `result$acf[s + 1, i, j]` uses
#' component \eqn{i} at time \eqn{t+s} and component \eqn{j} at time \eqn{t}.
#' The partial autoregression matrix in (16.5.44) is instead
#' \eqn{K_sV_s^{-1}}. Its entries need not lie in \eqn{[-1,1]}.
#' These definitions coincide for a univariate stationary process.
#'
#' \strong{Recursion.} With \eqn{C_h=\operatorname{Cov}(Z_{t+h},Z_t)}, initialize
#' \eqn{U_1=V_1=C_0}. Using forward and backward coefficients \eqn{F} and
#' \eqn{B} from the previous order, compute
#' \deqn{K_s=C_s-\sum_{k=1}^{s-1}F_{s-1,k}C_{s-k},\qquad
#' F_{s,s}=K_sV_s^{-1},\quad B_{s,s}=K_s'U_s^{-1}.}
#' For \eqn{1\leq k<s}, update both coefficient sets from the old order:
#' \deqn{F_{s,k}=F_{s-1,k}-F_{s,s}B_{s-1,s-k},\qquad
#' B_{s,k}=B_{s-1,k}-B_{s,s}F_{s-1,s-k}.}
#' The residual covariances for the next order are
#' \deqn{U_{s+1}=U_s-F_{s,s}K_s',\qquad
#' V_{s+1}=V_s-B_{s,s}K_s.}
#' This is the forward/backward recursion of (16.5.50)--(16.5.51), expressed
#' in R's covariance convention. The displayed inverses describe the mathematics;
#' implementations solve systems without explicitly forming matrix inverses.
#'
#' \strong{Estimation and numerical behavior.} Covariances are centered,
#' common-denominator estimates from `stats::acf()`, with denominator \eqn{n}.
#' Positive column rescaling improves numerical conditioning without changing
#' the estimator. For a valid joint covariance and positive residual variances,
#' individual coefficients are bounded by one in modulus, up to roundoff.
#' At lag one they equal ordinary lag-one cross-correlations.
#' Neither implementation clips correlations or adds a regularization term.
#' Singular or numerically singular regression covariances cause errors.
#' Short samples relative to the dimension and lag can therefore fail.
#' The last requested lag requires no further coefficient solve; lag one is
#' defined even for redundant columns having positive marginal variance.
#' Both implementations use \eqn{O(m^2 p)} storage and dense
#' \eqn{O(m^3 p^2)} arithmetic, where \eqn{p} is the maximum lag.
#' They can differ by floating point roundoff near singularity.
#'
#' This multivariate definition differs from the partial autoregression
#' coefficients returned by `stats::pacf()` and `stats::ar.yw(x)$partialacf`
#' in the R version for which this package was developed.
#' Loading this package changes neither function.
#'
#' @return An object of class `"acf"`, with the standard components `acf`, `lag`,
#'   `type`, `n.used`, `series` and `snames`. For a result named `result`,
#'   `result$acf` contains the coefficient array, with dimensions
#'   \eqn{p \times m \times m}, ordered by lag, component at time \eqn{t+s},
#'   and component at time \eqn{t}. Lags are in time
#'   units; plotting uses the standard signed cross-lag convention.
#'   Returned invisibly if `plot = TRUE`.
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
#' pc <- pacf_wei(x, lag.max = 10, plot = FALSE)  # C by default
#' pr <- pacf_wei(x, lag.max = 10, plot = FALSE, method = "R")
#' stopifnot(isTRUE(all.equal(pc$acf, pr$acf, tolerance = 1e-10)))
#' pc$acf
#' # Transpose both component axes to use Wei's theoretical convention:
#' aperm(pc$acf, c(1, 3, 2))
#' @export
#' @md
#' @useDynLib pacfwei, .registration = TRUE, .fixes = "C_"
pacf_wei <- function(x, lag.max = NULL, plot = TRUE,
                     na.action = stats::na.fail, ..., method = c("C", "R")) {
    method <- match.arg(method)
    series <- deparse1(substitute(x))
    x <- na.action(stats::as.ts(x))
    if (NROW(x) < 2L) stop("'x' must have at least two observations")
    x <- drop(x)
    if (!is.numeric(x) || is.complex(x)) stop("'x' must be real numeric")
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
