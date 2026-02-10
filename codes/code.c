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

int main()
{
	double h=0.01;
	int Nx=100;
	int Ny=100;
	
	double *phi=(double *) malloc(sizeof(double)*(Nx*Ny));
	double *psi=(double *) malloc(sizeof(double)*(Nx*Ny));
	double psi_top[Nx],psi_bottom[Nx];
	
	
	for(int i=0;i<Nx;i++)
	for(int j=0;j<Ny;j++)
	phi[i*Ny+j]=0.0;
	
	for(int i=0;i<Nx;i++)
	{
		psi_bottom[i]=10.0-10.*sin(2.*PI/Nx*i);
		psi_top[i]=-10.0+10.*sin(2.*PI/Nx*i);
	}
	
	
	
	poisson_2D_complex(h, Nx, Ny,psi_bottom,psi_top,phi,psi);
	FILE *F = fopen("../data/data2D.dat","w");
	
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Ny;j++)
		fprintf(F,"%f %f %f\n",i*h,j*h,psi[i+Nx*j]);
		fprintf(F,"\n");
	}
	
	fclose(F);
	return 0;
}
