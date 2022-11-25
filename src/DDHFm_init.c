#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void CentralDDHFT(double *sm, double *det, double *mu, double *sigma,
                int *nhalf, double *hft, double *factors, int *n, int *J);

extern void isotoneC(double *x, double *wt, int *nn, int *increasing);

static const R_CMethodDef CEntries[] = {
    {"CentralDDHFT", (DL_FUNC) &CentralDDHFT, 9},
    {"isotoneC",     (DL_FUNC) &isotoneC,     4},
    {NULL, NULL, 0}
};

void R_init_DDHFm(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

