/*
 * C port of Wei's partial lag correlation recursion, eqs. 16.5.43,50-51.
 * Based on Simone Giannerini's 2007 R implementation; GPL version 2 or later.
 *
 * The author of the original software is Simone Giannerini, Copyright (c) 2007.
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * https://www.R-project.org/Licenses/
 *
 * Input: REAL array (p+1,m,m), from stats::acf(type="covariance").
 * Output: REAL array (p,m,m), using R's later-endpoint/earlier-endpoint axes.
 * Covariance estimation, marginal rescaling, labels and plotting stay in R.
 * All dense products and solves use R's BLAS/LAPACK. No external dependencies.
 */
#define R_NO_REMAP
#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static R_xlen_t checked_product(R_xlen_t a, R_xlen_t b)
{
    if (a < 0 || b < 0 || (b && a > R_XLEN_T_MAX / b))
        Rf_error("covariance dimensions are too large");
    R_xlen_t n = a * b;
    if ((uintmax_t) n > SIZE_MAX / sizeof(double))
        Rf_error("workspace size exceeds addressable memory");
    return n;
}

static double *workspace(R_xlen_t n)
{
    return (double *) R_alloc((size_t) n, sizeof(double));
}

/* Z <- Z - A B.  Matrices are contiguous column-major m by m blocks. */
static void subtract_product(double *z, const double *a, const double *b, int m)
{
    const double minus = -1.0, one = 1.0;
    F77_CALL(dgemm)("N", "N", &m, &m, &m, &minus,
                   a, &m, b, &m, &one, z, &m FCONE FCONE);
}

/* out = B A^{-1}, by solving A' out' = B'. No explicit inverse is formed.
 * DGECON uses the LU factors to reject numerically singular systems, matching
 * R's solve.default tolerance of .Machine$double.eps in the reference code.
 */
static void right_solve(const double *b, const double *a, double *out,
                        int m, int lag, const char *side, double *lu,
                        double *rhs, int *pivot, double *work, int *iwork)
{
    double anorm = 0.0, rcond = 0.0;
    int info;
    for (int j = 0; j < m; ++j) {
        double sum = 0.0;
        for (int i = 0; i < m; ++i) {
            size_t ij = (size_t) i + (size_t) j * m;
            size_t ji = (size_t) j + (size_t) i * m;
            lu[ij] = a[ji];
            rhs[ij] = b[ji];
            sum += fabs(lu[ij]);
        }
        if (sum > anorm) anorm = sum;
    }
    F77_CALL(dgesv)(&m, &m, lu, &m, pivot, rhs, &m, &info);
    if (info > 0)
        Rf_error("singular %s residual covariance at lag %d (pivot %d)",
                 side, lag, info);
    if (info < 0)
        Rf_error("DGESV rejected argument %d at lag %d", -info, lag);
    F77_CALL(dgecon)("1", &m, lu, &m, &anorm, &rcond, work, iwork,
                    &info FCONE);
    if (info != 0)
        Rf_error("DGECON failed at lag %d (INFO=%d)", lag, info);
    if (!R_FINITE(rcond) || rcond < DBL_EPSILON)
        Rf_error("numerically singular %s residual covariance at lag %d (rcond=%g)",
                 side, lag, rcond);
    for (int j = 0; j < m; ++j)
        for (int i = 0; i < m; ++i)
            out[(size_t) i + (size_t) j*m] = rhs[(size_t) j + (size_t) i*m];
}

SEXP pacf_wei_cov(SEXP covariance)
{
    SEXP dim = Rf_getAttrib(covariance, R_DimSymbol);
    if (TYPEOF(covariance) != REALSXP || TYPEOF(dim) != INTSXP || XLENGTH(dim) != 3)
        Rf_error("'covariance' must be a double array with dimensions (p+1,m,m)");
    int nlags = INTEGER(dim)[0], m = INTEGER(dim)[1];
    if (nlags < 2 || m < 1 || INTEGER(dim)[2] != m)
        Rf_error("'covariance' must have dimensions (p+1,m,m), p>=1, m>=1");
    int p = nlags - 1;
    R_xlen_t mm = checked_product(m, m);
    R_xlen_t total = checked_product(mm, nlags);
    R_xlen_t nout = checked_product(mm, p);
    if (XLENGTH(covariance) != total)
        Rf_error("covariance length does not match its dimensions");
    const double *input = REAL(covariance);
    for (R_xlen_t i = 0; i < total; ++i) {
        if (!R_FINITE(input[i])) Rf_error("nonfinite covariance value");
        if ((i & 65535) == 0) R_CheckUserInterrupt();
    }
    /* Pack strided R lag slices into contiguous BLAS matrices. */
    double *g = workspace(total);
    for (int s = 0; s < nlags; ++s)
        for (R_xlen_t ij = 0; ij < mm; ++ij)
            g[(R_xlen_t)s*mm + ij] = input[s + (R_xlen_t)nlags*ij];

    SEXP outdim = PROTECT(Rf_allocVector(INTSXP, 3));
    INTEGER(outdim)[0] = p;
    INTEGER(outdim)[1] = INTEGER(outdim)[2] = m;
    SEXP answer = PROTECT(Rf_allocArray(REALSXP, outdim));
    double *output = REAL(answer);
    (void) nout; /* checked before allocation */
    double *vf = workspace(mm), *vb = workspace(mm), *cross = workspace(mm);
    memcpy(vf, g, (size_t)mm * sizeof(double));
    memcpy(vb, g, (size_t)mm * sizeof(double));

    /* At lag one no regression is needed, even if lag-zero covariance is
     * singular. All solver/history workspace can therefore be skipped. */
    R_xlen_t nhistory = p > 1 ? checked_product(mm, p-1) : 0;
    double *f = NULL, *b = NULL, *fnext = NULL, *bnext = NULL;
    double *fs = NULL, *bs = NULL, *lu = NULL, *rhs = NULL, *condwork = NULL;
    int *pivot = NULL, *iwork = NULL;
    if (p > 1) {
        f = workspace(nhistory); b = workspace(nhistory);
        fnext = workspace(nhistory); bnext = workspace(nhistory);
        fs = workspace(mm); bs = workspace(mm);
        lu = workspace(mm); rhs = workspace(mm);
        condwork = workspace(checked_product(4, m));
        pivot = (int *) R_alloc((size_t)m, sizeof(int));
        iwork = (int *) R_alloc((size_t)m, sizeof(int));
    }
    for (int s = 1; s <= p; ++s) {
        R_CheckUserInterrupt();
        memcpy(cross, g + (R_xlen_t)s*mm, (size_t)mm * sizeof(double));
        for (int k = 1; k < s; ++k) {
            if ((k & 63) == 0) R_CheckUserInterrupt();
            subtract_product(cross, f + (R_xlen_t)(k-1)*mm,
                             g + (R_xlen_t)(s-k)*mm, m);
        }
        for (int j = 0; j < m; ++j) {
            double db = vb[(size_t)j*(m+1)];
            if (!R_FINITE(db) || db <= 0.0)
                Rf_error("nonpositive backward residual variance at lag %d", s);
            for (int i = 0; i < m; ++i) {
                double df = vf[(size_t)i*(m+1)];
                if (!R_FINITE(df) || df <= 0.0)
                    Rf_error("nonpositive forward residual variance at lag %d", s);
                size_t ij = (size_t)i + (size_t)j*m;
                double value = cross[ij] / sqrt(df) / sqrt(db);
                if (!R_FINITE(value)) Rf_error("nonfinite partial correlation at lag %d", s);
                output[(s-1) + (R_xlen_t)p*ij] = value;
            }
        }
        if (s == p) break;
        right_solve(cross, vb, fs, m, s, "backward", lu, rhs, pivot, condwork, iwork);
        /* Second RHS is cross': reuse rhs as a temporary transposed block.
         * right_solve writes rhs, so bs initially holds the transposed input. */
        for (int j = 0; j < m; ++j)
            for (int i = 0; i < m; ++i)
                bs[(size_t)i+(size_t)j*m] = cross[(size_t)j+(size_t)i*m];
        right_solve(bs, vf, bs, m, s, "forward", lu, rhs, pivot, condwork, iwork);
        for (int k = 1; k < s; ++k) {
            if ((k & 63) == 0) R_CheckUserInterrupt();
            R_xlen_t offset = (R_xlen_t)(k-1)*mm;
            memcpy(fnext+offset, f+offset, (size_t)mm*sizeof(double));
            memcpy(bnext+offset, b+offset, (size_t)mm*sizeof(double));
            subtract_product(fnext+offset, fs, b+(R_xlen_t)(s-k-1)*mm, m);
            subtract_product(bnext+offset, bs, f+(R_xlen_t)(s-k-1)*mm, m);
        }
        memcpy(fnext+(R_xlen_t)(s-1)*mm, fs, (size_t)mm*sizeof(double));
        memcpy(bnext+(R_xlen_t)(s-1)*mm, bs, (size_t)mm*sizeof(double));
        double *swap = f; f = fnext; fnext = swap;
        swap = b; b = bnext; bnext = swap;
        const double minus = -1.0, one = 1.0;
        F77_CALL(dgemm)("N", "T", &m, &m, &m, &minus, fs, &m,
                       cross, &m, &one, vf, &m FCONE FCONE);
        subtract_product(vb, bs, cross, m);
        /* Symmetrize only roundoff: no clipping or regularization. */
        for (int j = 0; j < m; ++j)
            for (int i = 0; i < j; ++i) {
                size_t ij = (size_t)i+(size_t)j*m, ji = (size_t)j+(size_t)i*m;
                double v = 0.5*vf[ij]+0.5*vf[ji]; vf[ij]=vf[ji]=v;
                v = 0.5*vb[ij]+0.5*vb[ji]; vb[ij]=vb[ji]=v;
            }
    }
    UNPROTECT(2);
    return answer;
}
