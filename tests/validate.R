run_backend <- function(backend) {
# Usage: Rscript tests/validate.R [installed-library] [optional-Pinkham-rda]
# By default the Pinkham checks use the data bundled with the package.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 1L && nzchar(args[1L]))
    .libPaths(c(args[1L], .libPaths()))
stopifnot(requireNamespace("pacfwei", quietly = TRUE))
pacf_mts <- function(x, lag.max=NULL, convention=c("R","Wei")) {
  convention <- match.arg(convention)
  a <- pacfwei::pacf_wei(x,lag.max,plot=FALSE,method=backend)$acf
  if(convention=="Wei") aperm(a,c(1L,3L,2L)) else a
}
pacf_wei <- function(...) pacfwei::pacf_wei(...,method=backend)
native.symbol <- get("C_pacf_wei", asNamespace("pacfwei"))
native <- function(cov) .Call(native.symbol, cov)
prepare <- getFromNamespace(".pacf_prepare", "pacfwei")
original.stats.default <- getS3method("pacf", "default")
nchecks <- 0L
check <- function(ok, label) {
    if (!isTRUE(ok)) stop("FAILED: ", label, call. = FALSE)
    nchecks <<- nchecks + 1L
    cat("PASS ", label, "\n", sep = "")
}
near <- function(a, b, label, tolerance = 2e-12) {
    d <- if (length(a) == length(b)) max(abs(as.numeric(a)-as.numeric(b))) else Inf
    check(identical(dim(a), dim(b)) && is.finite(d) && d <= tolerance,
          sprintf("%s (max error %.3g)", label, d))
    invisible(d)
}
throws <- function(expr, label) {
    check(inherits(tryCatch({force(expr); NULL}, error = identity), "error"), label)
}

# Definition-based oracle: project the endpoints off ALL intervening vector
# coordinates using a block Toeplitz covariance's Schur complement. This
# deliberately does not use a Levinson recursion or the native implementation.
schur_reference <- function(x, p) {
    x <- as.matrix(x)
    m <- ncol(x)
    cc <- stats::acf(x, type="covariance", plot=FALSE, lag.max=p)$acf
    block <- function(i,j) {
        g <- matrix(cc[abs(i-j)+1L,,],m,m)
        if (i < j) t(g) else g
    }
    out <- array(NA_real_,c(p,m,m))
    for (s in seq_len(p)) {
        H <- do.call(rbind,lapply(0:s,function(i)
            do.call(cbind,lapply(0:s,function(j) block(i,j)))))
        endpoints <- c(seq_len(m),s*m+seq_len(m))
        S <- H[endpoints,endpoints,drop=FALSE]
        if (s > 1L) {
            ii <- (m+1L):(s*m)
            S <- S-H[endpoints,ii,drop=FALSE] %*%
                solve(H[ii,ii,drop=FALSE],H[ii,endpoints,drop=FALSE])
        }
        early <- seq_len(m)
        late <- m+early
        out[s,,] <- S[late,early,drop=FALSE] /
            outer(sqrt(diag(S)[late]),sqrt(diag(S)[early]))
    }
    out
}

cat(R.version.string,"\n")
set.seed(10)
scale.example <- cbind(x=rnorm(1000,sd=10000),y=rnorm(1000))
new <- pacf_mts(scale.example,10)
old <- stats::pacf(scale.example,10,plot=FALSE)$acf
cat(sprintf("stats::pacf range [%.12f, %.12f]; C Wei range [%.12f, %.12f]\n",
            min(old),max(old),min(new),max(new)))
check(max(abs(old)) > 500,"reported stats::pacf behavior reproduced")
check(max(abs(new)) <= 1,"scale example gives proper partial correlations")
near(new,schur_reference(scale.example,10),"scale example vs Schur projection")

set.seed(137)
x <- matrix(rnorm(1500),500,3,dimnames=list(NULL,c("alpha","beta","gamma")))
x[,2] <- .8*x[,1]+x[,2]
x[2:500,3] <- .5*x[1:499,1]+x[2:500,3]
p <- 12L
value <- pacf_mts(x,p)
near(value,schur_reference(x,p),"coupled three-series data vs Schur projection")
near(value[1,,],stats::acf(x,1,plot=FALSE)$acf[2,,],"lag one equals ordinary correlation")
near(pacf_mts(x,p,convention="Wei"),aperm(value,c(1L,3L,2L)),"Wei convention transposes lag matrices")
near(pacf_mts(x[nrow(x):1L,],p),aperm(value,c(1L,3L,2L)),"time reversal transposes lag matrices")
near(pacf_mts(sweep(x,2,c(.003,1e4,9),"*"),p),value,"positive rescaling invariance")
near(pacf_mts(sweep(x,2,c(1e-120,1e120,2),"*"),p),value,"extreme finite rescaling invariance")
scales <- c(-7,300,-.02)
signed <- value
for (s in seq_len(p)) signed[s,,] <- signed[s,,]*outer(sign(scales),sign(scales))
near(pacf_mts(sweep(x,2,scales,"*"),p),signed,"negative scales change row and column signs")
permutation <- c(3L,1L,2L)
near(pacf_mts(x[,permutation],p),value[,permutation,permutation],"series permutation equivariance")
for (j in seq_len(ncol(x)))
    near(drop(pacf_mts(x[,j,drop=FALSE],p)),
         drop(stats::pacf(x[,j],p,plot=FALSE)$acf),paste("univariate agreement, series",j))
near(pacf_mts(as.data.frame(x),p),value,"numeric data-frame input")
near(pacf_mts(x,1),value[1,,,drop=FALSE],"lag.max one")
near(pacf_mts(x,3.9),value[1:3,,,drop=FALSE],"fractional lag.max floors")
check(identical(dim(pacf_mts(x[1:15,1:2],100)),c(14L,2L,2L)),"lag.max caps at n-1")
check(identical(dimnames(value),list(as.character(seq_len(p)),colnames(x),colnames(x))),
      "lag and both series axes are named")
check(max(abs(stats::ar.yw(scale(x),order.max=p)$partialacf-value)) > .01,
      "standardizing before ar.yw is not the Wei estimator")

for (bad in list(0,-1,NA_real_,Inf,numeric(),c(1,2),"2"))
    throws(pacf_mts(x,bad),paste0("invalid lag.max rejected: ",paste(bad,collapse=",")))
for (bad in c(NA_real_,NaN,Inf,-Inf)) {
    xx <- x
    xx[3,2] <- bad
    throws(pacf_mts(xx,2),paste("nonfinite input rejected:",bad))
}
throws(pacf_mts(cbind(x[,1],1),2),"constant component rejected")
throws(pacf_mts(cbind(x[,1],0),2),"zero component rejected")
throws(pacf_mts(matrix(1,1,2),1),"fewer than two observations rejected")
throws(pacf_mts(matrix(numeric(),5,0),1),"zero series count rejected")
throws(pacf_mts(matrix(letters[1:10],5,2),2),"character input rejected")
throws(pacf_mts(x+1i,2),"complex input rejected")
duplicate <- cbind(x[,1],x[,1])
near(pacf_mts(duplicate,1)[1,,],stats::acf(duplicate,1,plot=FALSE)$acf[2,,],
     "redundant series permitted at lag one")
throws(pacf_mts(duplicate,2),"singular regression covariance fails explicitly")

cc <- prepare(x,3)$cov
near(native(cc),pacf_mts(x,3),"registered native routine matches R API")
throws(native(as.integer(cc)),"native rejects integer instead of real covariance")
throws(native(as.numeric(cc)),"native rejects missing dimensions")
throws(native(matrix(0,3,3)),"native rejects matrix instead of three-dimensional array")
throws(native(array(0,c(3,2,3))),"native rejects nonsquare covariance slices")
throws(native(array(0,c(1,2,2))),"native rejects covariance sequence without a positive lag")
throws(native(array(numeric(),c(3,0,0))),"native rejects zero-dimensional series")
cc.bad <- cc
cc.bad[1] <- NA_real_
throws(native(cc.bad),"native rejects NA covariance")
cc.bad[1] <- Inf
throws(native(cc.bad),"native rejects infinite covariance")
cc.bad <- cc
cc.bad[1,1,1] <- -1
throws(native(cc.bad),"native rejects nonpositive marginal variance")
ill.conditioned <- array(0,c(3L,2L,2L))
ill.conditioned[1,,] <- diag(c(1,1e-20))
condition.error <- tryCatch({native(ill.conditioned); NULL},error=identity)
check(inherits(condition.error,"error") &&
      grepl("numerically singular",conditionMessage(condition.error)),
      "native rejects invertible but numerically singular regression covariance")

xt <- stats::ts(x,start=c(2000,2),frequency=12)
out <- pacf_wei(xt,p,plot=FALSE)
check(inherits(out,"acf"),"wrapper returns acf class")
check(identical(out$snames,colnames(x)) && identical(out$n.used,nrow(x)),
      "wrapper retains series names and observation count")
near(out$acf,value,"wrapper matches array API")
signs <- matrix(1,3,3)
signs[lower.tri(signs)] <- -1
near(out$lag,outer(seq_len(p),signs/12),"wrapper retains frequency and signed cross-lags")
near(pacf_wei(x[,1,drop=FALSE],p,plot=FALSE)$acf,
     stats::pacf(x[,1],p,plot=FALSE)$acf,"one-column wrapper preserves univariate path")
xx <- xt
xx[1,] <- NA
throws(pacf_wei(xx,2,plot=FALSE),"wrapper rejects missing rows by default")
near(pacf_wei(xx,2,plot=FALSE,na.action=stats::na.omit)$acf,pacf_mts(xx[-1,],2),
     "wrapper applies NA omission action")
png.path <- tempfile("pacfwei-plot-",fileext=".png")
grDevices::png(png.path,width=900,height=900)
visible <- withVisible(pacf_wei(xt,3,plot=TRUE))
invisible(grDevices::dev.off())
check(!visible$visible && inherits(visible$value,"acf"),"plot path returns acf object invisibly")
unlink(png.path)
check(identical(getS3method("pacf","default"),original.stats.default),
      "stats::pacf remains unchanged")

data.env <- new.env()
utils::data("pinkham", package = "pacfwei", envir = data.env)
check(is.data.frame(data.env$pinkham) &&
      identical(dim(data.env$pinkham), c(54L, 2L)) &&
      identical(names(data.env$pinkham), c("advertising", "sales")) &&
      all(vapply(data.env$pinkham, is.integer, logical(1))) &&
      !anyNA(data.env$pinkham), "bundled Pinkham data have the documented format")
check(identical(as.integer(unlist(data.env$pinkham[1, ])), c(608L, 1016L)) &&
      identical(as.integer(unlist(data.env$pinkham[54, ])), c(564L, 1289L)),
      "bundled Pinkham endpoints are preserved")
if (length(args) >= 2L && file.exists(args[2L])) {
    supplied <- new.env()
    load(args[2L], envir = supplied)
    check(identical(data.env$pinkham, supplied$pinkham),
          "bundled Pinkham data are identical to supplied data")
}
{
    pinkham.diff <- apply(data.env$pinkham,2,diff)
    pinkham <- pacf_mts(pinkham.diff,6)
    wei <- pacf_mts(pinkham.diff,6,convention="Wei")
    near(pinkham,schur_reference(pinkham.diff,6),"Pinkham data vs Schur projection")
    near(round(wei[1,,],2),matrix(c(.05,.31,.22,.43),2,2),
         "Pinkham lag one, equation-defined Wei orientation",0)
    near(round(wei[2,,],2),matrix(c(-.53,-.32,-.28,-.10),2,2),
         "Pinkham lag two, equation-defined Wei orientation",0)
    cat("Table 16.2 has an apparent orientation inconsistency: its lag one\n",
        "matches R orientation, while its lags two through six match Wei orientation.\n",sep="")
}
cat(sprintf("\nAll %d validation checks passed.\n",nchecks))
}
for (backend in c("C","R")) {cat("\nBACKEND: ",backend,"\n");run_backend(backend)}
stopifnot(identical(getNamespaceExports("pacfwei"),"pacf_wei"))
set.seed(42); z <- matrix(rnorm(600),200,3)
stopifnot(identical(pacfwei::pacf_wei(z,5,plot=FALSE)$acf,
                    pacfwei::pacf_wei(z,5,plot=FALSE,method="C")$acf))
stopifnot(inherits(tryCatch(pacfwei::pacf_wei(z,method="bad"),error=identity),"error"))
for(method in c("C","R")) {
  stopifnot(inherits(tryCatch(pacfwei::pacf_wei(matrix(1:3,1,3),plot=FALSE,method=method),error=identity),"error"))
  set.seed(111); short <- matrix(rnorm(100*50),100,50)
  stopifnot(inherits(tryCatch(pacfwei::pacf_wei(short,50,plot=FALSE,method=method),error=identity),"error"))
}
cat("Single export, default C selection, invalid method and singular-backend checks passed.\n")
