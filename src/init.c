#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

extern SEXP pacf_wei_cov(SEXP);

static const R_CallMethodDef call_methods[] = {
    {"pacf_wei", (DL_FUNC) &pacf_wei_cov, 1},
    {NULL, NULL, 0}
};

void attribute_visible R_init_pacfwei(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, call_methods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

