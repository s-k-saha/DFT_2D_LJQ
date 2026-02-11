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




/*
void poisson_2D_complex(double h, int Nx, int Ny,double *psi_bottom, double *psi_top,double *phi, double *psi)
{
    int K = Nx/2 + 1;

    double *phi_1D = fftw_alloc_real(Nx);
    fftw_complex *phi_k_1D = fftw_alloc_complex(K);
    fftw_complex *phi_k = fftw_alloc_complex(K * Ny);

    fftw_plan plan_r2c = fftw_plan_dft_r2c_1d(Nx, phi_1D, phi_k_1D, FFTW_ESTIMATE);

    for(int j=0; j<Ny; j++)
    {
        // Copy row j into 1D buffer
        for(int i=0; i<Nx; i++)
            phi_1D[i] = phi[i + Nx*j];   // ✅ correct indexing

        fftw_execute(plan_r2c);

        // Store result
        for(int i=0; i<K; i++) {
            phi_k[j*K + i][0] = phi_k_1D[i][0];
            phi_k[j*K + i][1] = phi_k_1D[i][1];
        }
    }

    fftw_destroy_plan(plan_r2c);
    fftw_free(phi_1D);
    fftw_free(phi_k_1D);
		
		
		fftw_complex *psi_k_top = fftw_alloc_complex(K);
		fftw_complex *psi_k_bottom = fftw_alloc_complex(K);
		fftw_plan plan_r2c1= fftw_plan_dft_r2c_1d(Nx, psi_top, psi_k_top, FFTW_ESTIMATE);		
		fftw_plan plan_r2c2= fftw_plan_dft_r2c_1d(Nx, psi_bottom, psi_k_bottom, FFTW_ESTIMATE);
		
		fftw_execute(plan_r2c1);
		fftw_execute(plan_r2c2);
		fftw_destroy_plan(plan_r2c1);
		fftw_destroy_plan(plan_r2c2);
		
		double complex phi_1Dy[Ny];
		double complex psi_1Dy[Ny];
		
		fftw_complex *psi_k = fftw_alloc_complex(K * Ny);
		
		double k=0.;
		double lambda_k=0.;
		for(int i=0;i<K;i++)
		{
			k=2.*PI*i/(h*Nx);
			lambda_k=-k*k;
			
			for(int j=0;j<Ny;j++)
			phi_1Dy[j]=phi_k[j*K+i];
			
			
			poisson_1D_complex( h, Ny,psi_k_bottom[i],psi_k_top[i],phi_1Dy,psi_1Dy,lambda_k);
			
			for(int j=0;j<Ny;j++)
			psi_k[i+K*j]=psi_1Dy[j];
		}
		
		fftw_complex *psi_k_1D = fftw_alloc_complex(K);
		double *psi_1D = fftw_alloc_real(Nx);
		
		fftw_plan plan_c2r = fftw_plan_dft_c2r_1d(Nx, psi_k_1D, psi_1D, FFTW_ESTIMATE);
		
		for(int j=0;j<Ny;j++)
		{
			
				for(int i=0;i<K;i++)
				psi_k_1D[i] = psi_k[i+K*j];
				
				fftw_execute(plan_c2r);
				
				for(int i=0;i<Nx;i++)
				psi[i+j*Nx]=psi_1D[i]/Nx;
		}
    
    
    fftw_destroy_plan(plan_c2r);
    fftw_free(psi_k_1D);
    
    fftw_free(phi_k);
    fftw_free(psi_k);
    
    
    
    
}

*/

void poisson_2D_complex2(double h, int Nx, int Ny,
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

        double kx = 2.0 * M_PI * i / (Nx * h);
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




void main()
{
	double h=0.01;
	int Nx=100;
	int Ny=100;
	
	double phi[Nx*Ny];
	double psi[Nx*Ny];
	double psi_top[Nx],psi_bottom[Nx];
	
	
	for(int i=0;i<Nx;i++)
	for(int j=0;j<Ny;j++)
	phi[i*Ny+j]=0.0;
	
	for(int i=0;i<Nx;i++)
	{
		psi_bottom[i]=10.0;
		psi_top[i]=0.0;
	}
	
	
	
	poisson_2D_complex2(h, Nx, Ny,psi_bottom,psi_top,phi,psi);
	
	/*
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
	
	*/
	
	FILE *F = fopen("data2D.dat","w");
	
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		fprintf(F,"%f %f %f\n",i*h,j*h,psi[i+Nx*j]);
		fprintf(F,"\n");
	}
	
	fclose(F);
}

