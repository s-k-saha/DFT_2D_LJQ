#ifndef UTILITY_H
#define UTILITY_H

#include <complex.h>

/* -------------------------------------------------
   Constants
   ------------------------------------------------- */
#ifndef PI
#define PI 3.14159265358979
#endif

/* -------------------------------------------------
   Function prototypes
   ------------------------------------------------- */

/* 1D Poisson solver (complex-valued, Dirichlet BCs) */
void poisson_1D_complex(double h, int N,
                        double complex psi_left,
                        double complex psi_right,
                        double complex *phi,
                        double complex *psi,
                        double lambda_k);

/* 2D Poisson solver using FFT in x and FD in y */
void poisson_2D_complex(double h, int Nx, int Ny,
                        double *psi_bottom,
                        double *psi_top,
                        double *phi,
                        double *psi);

#endif /* UTILITY_H */

