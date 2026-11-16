/*
 * C port of Wei's partial autocorrelation recursion, eqs. 16.5.43,50-51.
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
 * Input: REAL covariance array (p+1,m,m) and REAL matrix A, with X = W A.
 * Covariances are estimated directly from whitened observations W. Each
 * residual covariance is mapped back as A' M A before normalization.
 * Output: REAL array (p,m,m), using R's component-at-t+s/component-at-t axes.
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
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>

static R_xlen_t product(R_xlen_t a, R_xlen_t b)
{
    if (a < 0 || b < 0 || (b && a > R_XLEN_T_MAX/b))
        Rf_error("covariance dimensions are too large");
    R_xlen_t n=a*b;
    if ((uintmax_t)n > SIZE_MAX/sizeof(double))
        Rf_error("workspace size exceeds addressable memory");
    return n;
}
static double *workspace(R_xlen_t n) { return (double *)R_alloc((size_t)n,sizeof(double)); }
/* Same LAPACK one-norm used by base::norm(A, "1"). DLANGE does not
 * reference WORK for this norm, but a valid address is supplied. */
static double norm1(const double *a, int m)
{
    double unused;
    return F77_CALL(dlange)("1", &m, &m, a, &m, &unused FCONE);
}
static void symmetrize(double *a,int m)
{
    for(int j=0;j<m;++j) for(int i=0;i<j;++i) {
        size_t ij=(size_t)i+(size_t)j*m,ji=(size_t)j+(size_t)i*m;
        double v=0.5*a[ij]+0.5*a[ji]; a[ij]=a[ji]=v;
    }
}
/* Z <- Z - A B. */
static void subtract_product(double *z,const double *a,const double *b,int m)
{
    const double minus=-1,one=1;
    F77_CALL(dgemm)("N","N",&m,&m,&m,&minus,a,&m,b,&m,&one,z,&m FCONE FCONE);
}
/* Original observations X = W A; residual covariances transform as A' M A. */
static void map_cov(const double *x,const double *a,double *out,double *tmp,int m)
{
    const double one=1,zero=0;
    F77_CALL(dgemm)("N","N",&m,&m,&m,&one,x,&m,a,&m,&zero,tmp,&m FCONE FCONE);
    F77_CALL(dgemm)("T","N",&m,&m,&m,&one,a,&m,tmp,&m,&zero,out,&m FCONE FCONE);
}
/* PSD is sufficient at exact boundaries, including the final requested lag. */
static void check_joint(const double *u,const double *v,const double *k,int m,int lag,
                        double *joint,double *eig,double *work,int lwork)
{
    int d=2*m,info;
    for(int j=0;j<m;++j) for(int i=0;i<m;++i) {
        size_t ij=(size_t)i+(size_t)j*m;
        joint[(size_t)i+(size_t)j*d]=u[ij];
        joint[(size_t)(i+m)+(size_t)(j+m)*d]=v[ij];
        joint[(size_t)i+(size_t)(j+m)*d]=k[ij];
        joint[(size_t)(j+m)+(size_t)i*d]=k[ij];
    }
    for(R_xlen_t i=0;i<product(d,d);++i)
        if(!R_FINITE(joint[i])) Rf_error("nonfinite residual covariance at lag %d",lag);
    double norm=norm1(joint,d);
    if(!R_FINITE(norm)) Rf_error("nonfinite residual covariance norm at lag %d",lag);
    double tol=64.0*d*DBL_EPSILON*fmax(1.0,norm);
    F77_CALL(dsyev)("N","U",&d,joint,&d,eig,work,&lwork,&info FCONE FCONE);
    if(info) Rf_error("residual covariance eigenvalue calculation failed at lag %d",lag);
    if(!R_FINITE(eig[0]) || eig[0]<-tol)
        Rf_error("invalid joint residual covariance at lag %d",lag);
}
/* Out = B A^-1, by Cholesky solves of A Out' = B'. No explicit inverse. */
static void right_solve(const double *b,const double *a,double *out,int m,int lag,
                        const char *side,double *factor,double *rhs,double *work,int *iwork)
{
    int info; double norm=norm1(a,m),rcond;
    memcpy(factor,a,(size_t)m*m*sizeof(double));
    F77_CALL(dpotrf)("U",&m,factor,&m,&info FCONE);
    if(info>0) Rf_error("nonpositive-definite %s residual covariance at lag %d",side,lag);
    if(info<0) Rf_error("DPOTRF rejected argument %d",-info);
    F77_CALL(dpocon)("U",&m,factor,&m,&norm,&rcond,work,iwork,&info FCONE);
    if(info) Rf_error("DPOCON failed at lag %d",lag);
    if(!R_FINITE(rcond) || rcond<sqrt(DBL_EPSILON))
        Rf_error("numerically singular %s residual covariance at lag %d",side,lag);
    for(int j=0;j<m;++j) for(int i=0;i<m;++i)
        rhs[(size_t)i+(size_t)j*m]=b[(size_t)j+(size_t)i*m];
    F77_CALL(dpotrs)("U",&m,&m,factor,&m,rhs,&m,&info FCONE);
    if(info) Rf_error("DPOTRS failed at lag %d",lag);
    for(int j=0;j<m;++j) for(int i=0;i<m;++i)
        out[(size_t)i+(size_t)j*m]=rhs[(size_t)j+(size_t)i*m];
}
SEXP pacf_wei_cov(SEXP covariance,SEXP map)
{
    SEXP dim=Rf_getAttrib(covariance,R_DimSymbol);
    if(TYPEOF(covariance)!=REALSXP || TYPEOF(dim)!=INTSXP || XLENGTH(dim)!=3)
        Rf_error("'covariance' must be a double array with dimensions (p+1,m,m)");
    int nlags=INTEGER(dim)[0],m=INTEGER(dim)[1];
    if(nlags<2 || m<1 || INTEGER(dim)[2]!=m || m>(INT_MAX-1)/6)
        Rf_error("invalid covariance dimensions");
    int p=nlags-1,d=2*m,lwork=3*d;
    R_xlen_t mm=product(m,m),total=product(mm,nlags),nout=product(mm,p);
    if(XLENGTH(covariance)!=total) Rf_error("covariance length does not match dimensions");
    SEXP md=Rf_getAttrib(map,R_DimSymbol);
    if(TYPEOF(map)!=REALSXP || TYPEOF(md)!=INTSXP || XLENGTH(md)!=2 ||
       INTEGER(md)[0]!=m || INTEGER(md)[1]!=m || XLENGTH(map)!=mm)
        Rf_error("'map' must be a finite double m by m matrix");
    const double *input=REAL(covariance),*a=REAL(map);
    for(R_xlen_t i=0;i<total;++i) {
        if(!R_FINITE(input[i])) Rf_error("nonfinite covariance value");
        if((i&65535)==0) R_CheckUserInterrupt();
    }
    for(R_xlen_t i=0;i<mm;++i) if(!R_FINITE(a[i])) Rf_error("'map' must contain finite values");
    double *g=workspace(total);
    for(int s=0;s<nlags;++s) for(R_xlen_t ij=0;ij<mm;++ij)
        g[(R_xlen_t)s*mm+ij]=input[s+(R_xlen_t)nlags*ij];
    double tol=64.0*m*DBL_EPSILON*fmax(1.0,norm1(g,m));
    for(int j=0;j<m;++j) for(int i=0;i<j;++i)
        if(fabs(g[(size_t)i+(size_t)j*m]-g[(size_t)j+(size_t)i*m])>tol)
            Rf_error("lag-zero covariance must be symmetric");
    symmetrize(g,m);
    SEXP outdim=PROTECT(Rf_allocVector(INTSXP,3));
    INTEGER(outdim)[0]=p;INTEGER(outdim)[1]=INTEGER(outdim)[2]=m;
    SEXP answer=PROTECT(Rf_allocArray(REALSXP,outdim));
    double *output=REAL(answer); (void)nout;
    double *vf=workspace(mm),*vb=workspace(mm),*cross=workspace(mm);
    memcpy(vf,g,(size_t)mm*sizeof(double));memcpy(vb,g,(size_t)mm*sizeof(double));
    double *mu=workspace(mm),*mv=workspace(mm),*mk=workspace(mm),*tmp=workspace(mm);
    double *joint=workspace(product(d,d)),*eig=workspace(d),*eigwork=workspace(lwork);
    R_xlen_t history=p>1?product(mm,p-1):0;
    double *f=NULL,*b=NULL,*fn=NULL,*bn=NULL,*fs=NULL,*bs=NULL;
    double *factor=NULL,*rhs=NULL,*condwork=NULL;int *iwork=NULL;
    if(p>1) {
        f=workspace(history);b=workspace(history);fn=workspace(history);bn=workspace(history);
        fs=workspace(mm);bs=workspace(mm);factor=workspace(mm);rhs=workspace(mm);
        condwork=workspace(product(3,m));iwork=(int *)R_alloc((size_t)m,sizeof(int));
    }
    double bound=64.0*m*DBL_EPSILON;
    for(int s=1;s<=p;++s) {
        R_CheckUserInterrupt();
        memcpy(cross,g+(R_xlen_t)s*mm,(size_t)mm*sizeof(double));
        for(int k=1;k<s;++k) {
            if((k&63)==0) R_CheckUserInterrupt();
            subtract_product(cross,f+(R_xlen_t)(k-1)*mm,g+(R_xlen_t)(s-k)*mm,m);
        }
        check_joint(vf,vb,cross,m,s,joint,eig,eigwork,lwork);
        map_cov(vf,a,mu,tmp,m);map_cov(vb,a,mv,tmp,m);map_cov(cross,a,mk,tmp,m);
        for(int j=0;j<m;++j) {
            double db=mv[(size_t)j*(m+1)];
            if(!R_FINITE(db)||db<=0) Rf_error("nonpositive backward residual variance at lag %d",s);
            for(int i=0;i<m;++i) {
                double df=mu[(size_t)i*(m+1)];
                if(!R_FINITE(df)||df<=0) Rf_error("nonpositive forward residual variance at lag %d",s);
                size_t ij=(size_t)i+(size_t)j*m;
                double value=mk[ij]/sqrt(df)/sqrt(db);
                if(!R_FINITE(value)) Rf_error("nonfinite partial correlation at lag %d",s);
                if(fabs(value)>1.0+bound) Rf_error("partial correlation exceeds numerical bounds at lag %d",s);
                output[s-1+(R_xlen_t)p*ij]=fmax(-1.0,fmin(1.0,value));
            }
        }
        if(s==p) break;
        right_solve(cross,vb,fs,m,s,"backward",factor,rhs,condwork,iwork);
        for(int j=0;j<m;++j) for(int i=0;i<m;++i)
            bs[(size_t)i+(size_t)j*m]=cross[(size_t)j+(size_t)i*m];
        right_solve(bs,vf,bs,m,s,"forward",factor,rhs,condwork,iwork);
        for(int k=1;k<s;++k) {
            if((k&63)==0) R_CheckUserInterrupt();
            R_xlen_t off=(R_xlen_t)(k-1)*mm;
            memcpy(fn+off,f+off,(size_t)mm*sizeof(double));memcpy(bn+off,b+off,(size_t)mm*sizeof(double));
            subtract_product(fn+off,fs,b+(R_xlen_t)(s-k-1)*mm,m);
            subtract_product(bn+off,bs,f+(R_xlen_t)(s-k-1)*mm,m);
        }
        memcpy(fn+(R_xlen_t)(s-1)*mm,fs,(size_t)mm*sizeof(double));
        memcpy(bn+(R_xlen_t)(s-1)*mm,bs,(size_t)mm*sizeof(double));
        double *swap=f;f=fn;fn=swap;swap=b;b=bn;bn=swap;
        const double minus=-1,one=1;
        F77_CALL(dgemm)("N","T",&m,&m,&m,&minus,fs,&m,cross,&m,&one,vf,&m FCONE FCONE);
        subtract_product(vb,bs,cross,m);symmetrize(vf,m);symmetrize(vb,m);
    }
    UNPROTECT(2);return answer;
}
