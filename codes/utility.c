#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include<limits.h>
#include <complex.h>
#include "utility.h"

#define PI 3.14159265358979


void poisson_1D_complex(double h, int N,
                        double complex psi_left,
                        double complex psi_right,
                        double complex *phi,
                        double complex *psi,
                        double lambda_k)
{
    int M = N - 2;  // number of interior unknowns

    if (M <= 0) {
        psi[0] = psi_left;
        psi[N-1] = psi_right;
        return;
    }

    /* Thomas algorithm coefficients (complex arithmetic) */
    double complex a =  1.0 / (h*h);          // subdiagonal
    double complex b = -2.0 / (h*h) + lambda_k; // diagonal
    double complex c =  1.0 / (h*h);          // superdiagonal

    /* Set Dirichlet boundary values */
    psi[0]   = psi_left;
    psi[N-1] = psi_right;

    /* Allocate Thomas arrays */
    double complex *bprime = (double complex*)malloc(M * sizeof(double complex));
    double complex *dprime = (double complex*)malloc(M * sizeof(double complex));

    /* ----------------------------
       Build RHS with boundary corrections
       ---------------------------- */
    dprime[0]   = phi[1] - a * psi_left;

    for (int i = 1; i < M-1; i++)
        dprime[i] = phi[i+1];

    dprime[M-1] = phi[N-2] - c * psi_right;

    /* ----------------------------
       Forward elimination
       ---------------------------- */
    bprime[0] = b;

    for (int i = 1; i < M; i++) {
        double complex m = a / bprime[i-1];
        bprime[i] = b - m * c;
        dprime[i] = dprime[i] - m * dprime[i-1];
    }

    /* ----------------------------
       Back substitution
       ---------------------------- */
    psi[N-2] = dprime[M-1] / bprime[M-1];

    for (int i = M-2; i >= 0; i--) {
        psi[i+1] = (dprime[i] - c * psi[i+2]) / bprime[i];
    }

    /* Cleanup */
    free(bprime);
    free(dprime);
}




void poisson_2D_complex(double h, int Nx, int Ny,
                        double *psi_bottom,
                        double *psi_top,
                        double *phi,
                        double *psi)
{
    int K = Nx/2 + 1;

    /* ---------- Forward FFT in x for RHS φ ---------- */

    double *phi_1D = fftw_alloc_real(Nx);
    fftw_complex *phi_k_1D = fftw_alloc_complex(K);
    fftw_complex *phi_k = fftw_alloc_complex(K * Ny);

    fftw_plan plan_r2c = fftw_plan_dft_r2c_1d(Nx, phi_1D, phi_k_1D, FFTW_ESTIMATE);

    for(int j = 0; j < Ny; j++) {
        for(int i = 0; i < Nx; i++)
            phi_1D[i] = phi[i + Nx*j];

        fftw_execute(plan_r2c);

        for(int i = 0; i < K; i++) {
            phi_k[j*K + i][0] = phi_k_1D[i][0];
            phi_k[j*K + i][1] = phi_k_1D[i][1];
        }
    }

    fftw_destroy_plan(plan_r2c);
    fftw_free(phi_1D);
    fftw_free(phi_k_1D);

    /* ---------- FFT boundary conditions in x ---------- */

    fftw_complex *psi_k_top = fftw_alloc_complex(K);
    fftw_complex *psi_k_bottom = fftw_alloc_complex(K);

    fftw_plan plan_bc_top = fftw_plan_dft_r2c_1d(Nx, psi_top, psi_k_top, FFTW_ESTIMATE);
    fftw_plan plan_bc_bot = fftw_plan_dft_r2c_1d(Nx, psi_bottom, psi_k_bottom, FFTW_ESTIMATE);

    fftw_execute(plan_bc_top);
    fftw_execute(plan_bc_bot);

    fftw_destroy_plan(plan_bc_top);
    fftw_destroy_plan(plan_bc_bot);

    /* ---------- Solve Poisson in y for each kx ---------- */

    fftw_complex *psi_k = fftw_alloc_complex(K * Ny);

    double complex *phi_1Dy = malloc(sizeof(double complex)*Ny);
    double complex *psi_1Dy = malloc(sizeof(double complex)*Ny);

    for(int i = 0; i < K; i++) {

        double kx = 2.0 * PI * i / (Nx * h);
        double lambda_k = -kx * kx;

        /* Extract RHS column */
        for(int j = 0; j < Ny; j++) {
            phi_1Dy[j] = phi_k[j*K + i][0] + I * phi_k[j*K + i][1];
        }

        double complex bc_bottom = psi_k_bottom[i][0] + I * psi_k_bottom[i][1];
        double complex bc_top    = psi_k_top[i][0]    + I * psi_k_top[i][1];

        poisson_1D_complex(h, Ny, bc_bottom, bc_top,
                           phi_1Dy, psi_1Dy, lambda_k);

        /* Store ψ̂(kx,y) */
        for(int j = 0; j < Ny; j++) {
            psi_k[j*K + i][0] = creal(psi_1Dy[j]);
            psi_k[j*K + i][1] = cimag(psi_1Dy[j]);
        }
    }

    free(phi_1Dy);
    free(psi_1Dy);
    fftw_free(psi_k_top);
    fftw_free(psi_k_bottom);
    fftw_free(phi_k);

    /* ---------- Inverse FFT in x ---------- */

    fftw_complex *psi_k_1D = fftw_alloc_complex(K);
    double *psi_1D = fftw_alloc_real(Nx);

    fftw_plan plan_c2r = fftw_plan_dft_c2r_1d(Nx, psi_k_1D, psi_1D, FFTW_ESTIMATE);

    for(int j = 0; j < Ny; j++) {

        for(int i = 0; i < K; i++) {
            psi_k_1D[i][0] = psi_k[j*K + i][0];
            psi_k_1D[i][1] = psi_k[j*K + i][1];
        }

        fftw_execute(plan_c2r);

        for(int i = 0; i < Nx; i++)
            psi[i + Nx*j] = psi_1D[i] / Nx;   // FFTW normalization
    }

    fftw_destroy_plan(plan_c2r);
    fftw_free(psi_k_1D);
    fftw_free(psi_1D);
    fftw_free(psi_k);
}





