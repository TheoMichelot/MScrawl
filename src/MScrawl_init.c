#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

// Generated with tools::package_native_routine_registration_skeleton

/* .Call calls */
extern SEXP _MScrawl_kalman_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _MScrawl_smooth_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_MScrawl_kalman_rcpp", (DL_FUNC) &_MScrawl_kalman_rcpp, 5},
    {"_MScrawl_smooth_rcpp", (DL_FUNC) &_MScrawl_smooth_rcpp, 5},
    {NULL, NULL, 0}
};

void R_init_MScrawl(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}