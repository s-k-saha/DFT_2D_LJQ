#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include<limits.h>
#include <complex.h>

#define PI 3.14159265358979
/* 
   fft_x_forward()  : FFT along x for each y
   fft_x_backward() : inverse FFT
   These are placeholders for FFTW or similar
*/

void poisson_1D(double h, int N,
                double psi_left, double psi_right,
                double *phi, double *psi,double lambda_k)
{
    /* number of interior unknowns */
    int M = N - 2;

    /* coefficients */
    double a =  1.0 / (h*h);   /* subdiagonal  */
    double b = -2.0 / (h*h) + lambda_k;   /* diagonal     */
    double c =  1.0 / (h*h);   /* superdiagonal*/

    /* set boundary values */
    psi[0]   = psi_left;
    psi[N-1] = psi_right;

    /* trivial case: no interior points */
    if (M <= 0)
        return;

    /* allocate Thomas arrays */
    double *bprime = (double *)malloc(M * sizeof(double));
    double *dprime = (double *)malloc(M * sizeof(double));

    /* ----------------------------
       Build RHS with BC corrections
       ---------------------------- */

    /* first interior point (i = 1) */
    dprime[0] = phi[1] - a * psi_left;

    /* middle points */
    for (int i = 1; i < M-1; i++)
        dprime[i] = phi[i+1];

    /* last interior point (i = N-2) */
    dprime[M-1] = phi[N-2] - c * psi_right;

    /* ----------------------------
       Thomas forward elimination
       ---------------------------- */

    bprime[0] = b;

    for (int i = 1; i < M; i++) {
        double m = a / bprime[i-1];
        bprime[i] = b - m * c;
        dprime[i] = dprime[i] - m * dprime[i-1];
    }

    /* ----------------------------
       Back substitution
       ---------------------------- */

    /* last unknown: psi[N-2] */
    psi[N-2] = dprime[M-1] / bprime[M-1];

    for (int i = M-2; i >= 0; i--) {
        psi[i+1] = (dprime[i] - c * psi[i+2]) / bprime[i];
    }

    /* cleanup */
    free(bprime);
    free(dprime);
}

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



void poisson_2D(double h, int Nx, int Ny,
                double *psi_bottom, double *psi_top,
                double phi[Nx][Ny], double psi[Nx][Ny])
{
    int K = Nx/2 + 1;       // Number of Fourier modes in x
    double Lx = h * (Nx - 1);

    // Flatten phi into 1D for FFTW
    double *phi2D = (double*)malloc(Nx*Ny*sizeof(double));
    double *psi2D = (double*)malloc(Nx*Ny*sizeof(double));

    for (int i=0;i<Nx;i++)
        for (int j=0;j<Ny;j++)
            phi2D[i*Ny + j] = phi[i][j];

    // FFT in x (real-to-complex)
    fftw_complex *phi_k = fftw_alloc_complex(K*Ny);
    fftw_plan plan_r2c = fftw_plan_dft_r2c_2d(Nx, Ny, phi2D, phi_k, FFTW_ESTIMATE);
    fftw_execute(plan_r2c);
    fftw_destroy_plan(plan_r2c);

    // FFT of boundary conditions in x
    fftw_complex *psi_bottom_k = fftw_alloc_complex(K);
    fftw_complex *psi_top_k    = fftw_alloc_complex(K);

    fftw_plan plan_b = fftw_plan_dft_r2c_1d(Nx, psi_bottom, psi_bottom_k, FFTW_ESTIMATE);
    fftw_execute(plan_b);
    fftw_destroy_plan(plan_b);

    fftw_plan plan_t = fftw_plan_dft_r2c_1d(Nx, psi_top, psi_top_k, FFTW_ESTIMATE);
    fftw_execute(plan_t);
    fftw_destroy_plan(plan_t);

    // Allocate storage for solution in Fourier space
    double complex *psi_k_1D = (double complex*)malloc(Ny*sizeof(double complex));
    double complex **psi_k = (double complex**)malloc(K * sizeof(double complex*));
    for (int i=0;i<K;i++)
        psi_k[i] = (double complex*)malloc(Ny*sizeof(double complex));

    // Solve 1D Poisson in y for each x-Fourier mode
    for (int i=0;i<K;i++)
    {
        double kx = 2.0*PI*i/Lx;      // wave number
        double lambda_k = -kx*kx;//-4./h/h*(sin(PI*i/(Nx-1))*sin(PI*i/(Nx-1)));     // Fourier transformed Laplacian in x

        // Extract phi_k along y for this mode
        for (int j=0;j<Ny;j++)
            psi_k_1D[j] = phi_k[i*Ny + j][0] + I*phi_k[i*Ny + j][1];

        // Extract complex boundary values for this mode
        double complex bc_bottom = psi_bottom_k[i][0] + I*psi_bottom_k[i][1];
        double complex bc_top    = psi_top_k[i][0] + I*psi_top_k[i][1];

        // Solve 1D complex Poisson along y
        poisson_1D_complex(h, Ny, bc_bottom, bc_top, psi_k_1D, psi_k[i], lambda_k);
    }

    // Flatten psi_k into fftw_complex array for inverse FFT
    fftw_complex *psi_k_flat = fftw_alloc_complex(K*Ny);
    for (int i=0;i<K;i++)
        for (int j=0;j<Ny;j++)
        {
            psi_k_flat[i*Ny + j][0] = creal(psi_k[i][j]);
            psi_k_flat[i*Ny + j][1] = cimag(psi_k[i][j]);
        }

    // Inverse FFT to get real-space solution
    fftw_plan plan_c2r = fftw_plan_dft_c2r_2d(Nx, Ny, psi_k_flat, psi2D, FFTW_ESTIMATE);
    fftw_execute(plan_c2r);
    fftw_destroy_plan(plan_c2r);

    // Copy solution back to 2D array and normalize by Nx (FFTW does not normalize)
    for (int i=0;i<Nx;i++)
        for (int j=0;j<Ny;j++)
            psi[i][j] = psi2D[i*Ny + j]/Nx;

    // Free memory
    fftw_free(phi_k);
    fftw_free(psi_bottom_k);
    fftw_free(psi_top_k);
    fftw_free(psi_k_flat);

    free(phi2D);
    free(psi2D);

    for (int i=0;i<K;i++)
        free(psi_k[i]);
    free(psi_k);

    free(psi_k_1D);
}


void main()
{
	double h=0.01;
	int Nx=100;
	int Ny=100;
	
	double phi[Nx][Ny];
	double psi[Nx][Ny];
	
	double psi_bottom[Nx];

	double psi_top[Nx];
	
	
	
	for(int i=0;i<Nx;i++)
	{psi_bottom[i]=10.0;psi_top[i]=0.0;}
	
	for(int i=0;i<Nx;i++)
	for(int j=0;j<Ny;j++)
	phi[i][j]=psi[i][j]=0.0;
	
	for(int i=0;i<Nx;i++)
	{psi[i][0]=psi_bottom[i];psi[i][Ny-1]=psi_top[i];}
	
	
	poisson_2D( h,  Nx, Ny,psi_bottom, psi_top, phi, psi);
	
	
	FILE *F = fopen("data2D.dat","w");
	
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		fprintf(F,"%f %f %f\n",i*h,j*h,psi[i]);
		fprintf(F,"\n");
	}
	
	fclose(F);
	
}

