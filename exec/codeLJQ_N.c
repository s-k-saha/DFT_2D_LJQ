#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <fftw3.h>
#include<limits.h>
#include <complex.h>
#include "utility.h"

#define MIN(x,y) ((x<y)? x : y)
#define MAX(x,y) ((x>y)? x : y)

#define PI 3.14159265358979
#define R 0.5
#define sigma 1
#define dz 0.025
#define dx 0.025
//#define mu .6045291013//-3.8214630890212753

#define rc 2.5
#define Nx 1999
#define Nz 1999

#define Lx (dx*(Nx-1))
#define Lz (dz*(Nz-1))

#define K (Nz/2+1)//Nx/2+1 
#define NiR (R/dz)
#define NiLJ ((int)(rc/dz))//rc/dz
#define R2 R*R
#define R3 R2*R

#define N (Nx*Nz) //Nx*Nz
#define NiR2 (4*NiR*NiR) //4*NiR*NiR

/*
full system: 0->solvent, 1-> +ion, 2-> -ion
invariant along x and y (rho(z))
full grand-canonical minimization
*/

double lambdaB;
double Vq;

const double alpha=0.01;
const double rmin=pow(2.,1./6);
//const double rc=2.5;
double eps=0.0;//0.9793868;
double eps_i=0.0;

double rhob_solvent_l=0.0;//0.6039361095188667
double rhob_solvent_g=0.0;//0.026202777905486873

double rhob_solute_l=0.0;//0.019845282259923474
double rhob_solute_g=0.0;//6.316260027217214e-05


double mu_solvent=0.0;
double mu_solute=0.0;

int Nbatch=0;

double ew;
double Nrho[3];

double rho[3][Nx*Nz],rhonew[3][Nx*Nz],rhocopy[Nx*Nz],Vext[Nx*Nz],Vext_Q[Nx*Nz];
double phi[Nx*Nz];
double psi[Nx*Nz];
double psi_bottom[Nx];
double psi_top[Nx];

double n0[Nx*Nz],n1[Nx*Nz],n2[Nx*Nz],n3[Nx*Nz],n1vx[Nx*Nz],n1vz[Nx*Nz],n2vx[Nx*Nz],n2vz[Nx*Nz];
double dphidn0[Nx*Nz],dphidn1[Nx*Nz],dphidn2[Nx*Nz],dphidn3[Nx*Nz],dphidn1vx[Nx*Nz],dphidn1vz[Nx*Nz],dphidn2vx[Nx*Nz],dphidn2vz[Nx*Nz];
double c1[3][Nx*Nz];
double c1_temp[Nx*Nz];

fftw_complex  omega3[Nx*K];
fftw_complex  omega2[Nx*K];
fftw_complex  omega1[Nx*K];
fftw_complex  omega0[Nx*K];
//double U[Nx][Nz];

fftw_complex omega1vx[Nx*K];
fftw_complex omega1vz[Nx*K];

fftw_complex omega2vx[Nx*K];
fftw_complex omega2vz[Nx*K];

fftw_complex  UfilterFFT[Nx*K];

// FFTW global buffers and plans
double *fft_in;
double *fft_in2; // only needed for conv_FFT2D
double *fft_out_real;

fftw_complex *fft_f;
fftw_complex *fft_g;
fftw_complex *fft_h;

fftw_plan plan_r2c_f;
fftw_plan plan_r2c_g;
fftw_plan plan_c2r;



double dist(int,int,int,int);

void getn(), getc1_fmt(),getc1_LJ(),filterc1(), getVext(), rhoinit(), iterate();
void getomega3(),getomega2(),getomega1(),getomega0(),getomega1v(),getomega2v(),rhocpy(),filterrho(),write_rho(double,int),conv_FFT2D_2(double*,fftw_complex *,double*),getc1_LJ2(),getc1_LJ2_i(int,int,double),getUFilterFFT(),read_params(),initialize_dataframes(int,double*),initialize_all_dataframes();

void init_fftw(),cleanup_fftw();

double aux_J5(double, double);
double aux_J11(double, double);
double J5(double, double,double);
double J11(double, double,double);
double ULJ(double);
void get_mu();

int getDx(int),getDz(int);
void add_c(double*,double*);
void add_cv(double*,double*);
void conv_FFT2D(double*, double*,double*);



double dist(int i1,int j1,int i2,int j2)
{
	double Dx=dx*abs(i2-i1);
	double Dz=dz*abs(j2-j1);
	Dx=MIN(Dx,Lx-Dx);
	Dz=MIN(Dz,Lz-Dz);
	return Dx*Dx+Dz*Dz;
}

double aux_J11(double r, double a)
{
	double t1=(315./1280) * (pow(a,-2))*(pow(r,-11));
	double t2=(210./1280) * (pow(a,-4))*(pow(r,-9));
	double t3=(168./1280) * (pow(a,-6))*(pow(r,-7));
	double t4=(144./1280) * (pow(a,-8))*(pow(r,-5));
	double t5=(1./10) * (pow(a,-10))*(pow(r,-3));
	
	double temp=0.;
	if(a==r)
		temp=0.;
	else
		temp=atan(sqrt(a*a/r/r -1));
	
	
	return (r*sqrt(a*a-r*r)*(t1+t2+t3+t4+t5)) + (315./1280)*temp*(pow(r,-11));
}

double aux_J5(double r, double a)
{
	double temp=0;
	if(a==r)
		temp=0.;
	else	
		temp=atan(sqrt(a*a/r/r -1));
	
	double t1=(3./8)*(pow(a,-2))*(pow(r,-5));
	double t2=(1./4)*(pow(a,-4))*(pow(r,-3));
		
	return (r*sqrt(a*a-r*r)*(t1+t2)) + (3./8)*(pow(r,-5))*temp;
	
}


double J11(double r, double a, double b)
{
	if(r>=1.0)
	return aux_J11(r,b)-aux_J11(r,a);
	
	
	double t1=(1./11) * (-pow(b,-11)+pow(a,-11));
	double t2=(1./26) * (-pow(b,-13)+pow(a,-13)) * pow(r,2);
	double t3= (1./40) * (-pow(b,-15)+pow(a,-15)) * pow(r,4);
	double t4=	(5./272) * (-pow(b,-17)+pow(a,-17)) * pow(r,6);
	double t5= (35./2432) * (-pow(b,-19)+pow(a,-19)) * pow(r,8);
	double t6= (3./256)* (-pow(b,-21)+pow(a,-21)) * pow(r,10);
	double t7=(231./23552)* (-pow(b,-23)+pow(a,-23)) * pow(r,12);
	double t8= (429./51200)* (-pow(b,-25)+pow(a,-25)) * pow(r,14);
  
  return (t1+t2+t3+t4+t5+t6+t7+t8);
}

double J5(double r, double a, double b)
{
	if(r>=1.0)
		return aux_J5(r,b)-aux_J5(r,a);
		
	double t1=(1./5) * (-pow(b,-5)+pow(a,-5));
	double t2=(1./14) * (-pow(b,-7)+pow(a,-7)) * pow(r,2);
	double t3= (1./24) * (-pow(b,-9)+pow(a,-9)) * pow(r,4);
	double t4=	(5./176) *  (-pow(b,-11)+pow(a,-11)) * pow(r,6);
	double t5= (35./1664) *  (-pow(b,-13)+pow(a,-13)) * pow(r,8);
	double t6= (21./1280)* (-pow(b,-15)+pow(a,-15)) * pow(r,10);
	double t7= (231./17408)* (-pow(b,-17)+pow(a,-17)) * pow(r,12);
	double t8= (429./38912)* (-pow(b,-19)+pow(a,-19)) * pow(r,14);
	
	return (t1+t2+t3+t4+t5+t6+t7+t8);
    
}



void getomega3()
{
	double kx=0.0,kz=0.0;
  double dkx=2.*PI/Lx,dkz=2.*PI/Lz;
  double k;
  for(int i=0;i<Nx;i++)
  {
  	if(i<=Nx/2)
  		kx=dkx*i;
  	else
  		kx=-dkx*(Nx-i);
  	for(int j=0;j<K;j++)
  	{
  		kz=dkz*j;
			k=sqrt(kx*kx+kz*kz);
			if(i==0 && j==0)
			{
				omega3[i*K+j][0]=4.*PI*R*R*R/3;
				omega3[i*K+j][1]=0.0;
			}
			else
			{
				omega3[i*K+j][0]=8./3*PI*R*R*(j1(R*k)/(k));
				omega3[i*K+j][1]=0.0;
  		}
  	}
  }
	
}

void getomega2()
{
  double kx=0.0,kz=0.0;
  double dkx=2.*PI/Lx,dkz=2.*PI/Lz;
  double k;
  for(int i=0;i<Nx;i++)
  {
  	if(i<=Nx/2)
  		kx=dkx*i;
  	else
  		kx=-dkx*(Nx-i);
  	for(int j=0;j<K;j++)
  	{
  		kz=dkz*j;
			k=sqrt(kx*kx+kz*kz);
			if(i==0 && j==0)
			{
				omega2[i*K+j][0]=4.*PI*R*R;
				omega2[i*K+j][1]=0.0;
			}
			else
			{
				omega2[i*K+j][0]=8.*PI*R*R*(j1(R*k)/(R*k));
				omega2[i*K+j][1]=0.0;
  		}
  	}
  }
}

void getomega1()
{
		for(int i=0;i<Nx;i++)
		for(int j=0;j<K;j++)
		{
		omega1[i*K+j][0]=omega2[i*K+j][0]*(1./4./PI/R);
		omega1[i*K+j][1]=omega2[i*K+j][1]*(1./4./PI/R);
		}
}
void getomega0()
{
		for(int i=0;i<Nx;i++)
		for(int j=0;j<K;j++)
		{
		omega0[i*K+j][0]=omega2[i*K+j][0]*(1./4./PI/R/R);
		omega0[i*K+j][1]=omega2[i*K+j][1]*(1./4./PI/R/R);
		}
}

void getomega2v()
{
	double kx=0.0,kz=0.0;
	double dkx=2.*PI/Lx,dkz=2.*PI/Lz;
	double k;
	for(int i=0;i<Nx;i++)
	{
		if(i<=Nx/2)
  		kx=dkx*i;
  	else
  		kx=-dkx*(Nx-i);
  	for(int j=0;j<K;j++)
  	{
  		kz=dkz*j;
			k=sqrt(kx*kx+kz*kz);
		
			omega2vx[i*K+j][1]=-kx*omega3[i*K+j][0];
			omega2vx[i*K+j][0]=0.0;
			omega2vz[i*K+j][1]=-kz*omega3[i*K+j][0];
			omega2vz[i*K+j][0]=0.0;	
		}
	
	}
}


void getomega1v()
{
	for(int i=0;i<Nx;i++)
	for(int j=0;j<K;j++)
	{
		omega1vx[i*K+j][0]=omega2vx[i*K+j][0]/(4.*PI*R);
		omega1vx[i*K+j][1]=omega2vx[i*K+j][1]/(4.*PI*R);
		
		omega1vz[i*K+j][0]=omega2vz[i*K+j][0]/(4.*PI*R);
		omega1vz[i*K+j][1]=omega2vz[i*K+j][1]/(4.*PI*R);
		
	}
}


void get_mu()
{
	double rhob=rhob_solvent_l+2.*rhob_solute_l;
	double eta = (rhob)*(4./3*PI*R3);
	mu_solvent=log(rhob_solvent_l) + (14.*eta-13.*eta*eta+5.*eta*eta*eta)/(2.*(1.-eta)*(1.-eta)*(1.-eta)) - log(1.-eta) - (1.171861897)*(4.*PI*rhob_solvent_l)*eps - (1.171861897)*(8.*PI*rhob_solute_l)*eps_i;
	
	mu_solute=log(rhob_solute_l) + (14.*eta-13.*eta*eta+5.*eta*eta*eta)/(2.*(1.-eta)*(1.-eta)*(1.-eta)) - log(1.-eta) - (1.171861897)*
	(4.*PI*rhob_solvent_l)*eps_i;
	printf("mu_solvent:%f\nmu_solute:%f\n",mu_solvent,mu_solute);
}

void rhoinit()
{
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			{
				rho[0][i*Nz+j]=rhob_solvent_g;
				rho[1][i*Nz+j]=rhob_solute_g;
				rho[2][i*Nz+j]=rhob_solute_g;
			}
	for(int i=0;i<Nx;i++) //Setting walls left 
		for(int j=0;j<NiR;j++)
			rho[0][i*Nz+j]=rho[1][i*Nz+j]=rho[2][i*Nz+j]=0.0;
	
	int xbuf=Nx*.25;
	int zbuf=Nz*.4;
	for(int i=xbuf;i<Nx-xbuf;i++)
		for(int j=NiR;j<NiR+zbuf;j++)
			{
				rho[0][i*Nz+j]=rhob_solvent_l;
				rho[1][i*Nz+j]=rhob_solute_l;
				rho[2][i*Nz+j]=rhob_solute_l;
			}
	Nrho[0]=Nrho[1]=Nrho[2]=0;
	for(int i=0;i<Nx;i++)
	for(int j=NiR;j<Nz;j++)
		{
			Nrho[0]+=(rho[0][i*Nz+j]-rhob_solvent_g);
			Nrho[1]+=(rho[1][i*Nz+j]-rhob_solute_g);
			Nrho[2]+=(rho[2][i*Nz+j]-rhob_solute_g);
			
		}
	
	for(int j=0;j<3;j++)
		Nrho[j]=(((Nrho[j]*dx*dz)/Lx));

	printf("N_solvent: %f\nN_solute: %f\n",Nrho[0],Nrho[1]);
	/*
	for(int i=0;i<Nx;i++)
		for(int j=Nz-NiR+1;j<Nz;j++) //Setting walls right
			rho[i*Nz+j]=0.0;
	*/
}


void getn()
{

	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			{
				n0[i1*Nz+j1]=0;
				n1[i1*Nz+j1]=0;
				n2[i1*Nz+j1]=0;
				n3[i1*Nz+j1]=0;
				n1vx[i1*Nz+j1]=0;
				n1vz[i1*Nz+j1]=0;
				n2vx[i1*Nz+j1]=0;
				n2vz[i1*Nz+j1]=0;
			}
	
	conv_FFT2D_2(rhocopy, omega0, n0);
	conv_FFT2D_2(rhocopy, omega1, n1);
	conv_FFT2D_2(rhocopy, omega2,n2);
	conv_FFT2D_2(rhocopy, omega3, n3);
	conv_FFT2D_2(rhocopy, omega1vx, n1vx);
	conv_FFT2D_2(rhocopy, omega1vz, n1vz);
	conv_FFT2D_2(rhocopy, omega2vx, n2vx);
	conv_FFT2D_2(rhocopy, omega2vz, n2vz);
	
}

void read_params()
{
	FILE *fp=fopen("params.txt", "r");
	char line[256];
	while (fgets(line, sizeof(line), fp)) 
	{
    if (sscanf(line, "lambdaB=%lf", &lambdaB) == 1) continue;
    if (sscanf(line, "Vq=%lf", &Vq) == 1) continue;
    if (sscanf(line, "ew=%lf", &ew) == 1) continue;
    if (sscanf(line, "eps=%lf", &eps) == 1) continue;
    if (sscanf(line, "eps_i=%lf", &eps_i) == 1) continue;
    if (sscanf(line, "rhob_solvent_l=%lf", &rhob_solvent_l) == 1) continue;
    if (sscanf(line, "rhob_solvent_g=%lf", &rhob_solvent_g) == 1) continue;
    if (sscanf(line, "rhob_solute_l=%lf", &rhob_solute_l) == 1) continue;
    if (sscanf(line, "rhob_solute_g=%lf", &rhob_solute_g) == 1) continue; 
    if (sscanf(line, "Nbatch=%d", &Nbatch) == 1) continue;
  }
	
  printf("\n----------------\nparams:\nNbatch:%d\nlambdaB:%f\nVq:%f\new:%f\neps:%f\neps_i:%f\nrhob_solvent_l:%f\nrhob_solvent_g:%f\nrhob_solute_l:%f\nrhob_solute_g:%f\n----------------\n",Nbatch,lambdaB,Vq,ew,eps,eps_i,rhob_solvent_l,rhob_solvent_g,rhob_solute_l,rhob_solute_g); 
  fclose(fp);
}

void conv_FFT2D(double *f, double *g, double *h2)
{
    memcpy(fft_in,  f, sizeof(double) * N);
    memcpy(fft_in2, g, sizeof(double) * N);

    fftw_execute(plan_r2c_f);
    fftw_execute(plan_r2c_g);

    for (int i = 0; i < Nx * K; ++i) {
        double a = fft_f[i][0], b = fft_f[i][1];
        double c = fft_g[i][0], d = fft_g[i][1];

        fft_h[i][0] = a*c - b*d;
        fft_h[i][1] = a*d + b*c;
    }

    fftw_execute(plan_c2r);

    for (int i = 0; i < N; i++)
        h2[i] = (fft_out_real[i] / N) * dx * dz;
}

void conv_FFT2D_2(double *f, fftw_complex *fft_g_input, double *h2)
{
    memcpy(fft_in, f, sizeof(double) * N);

    fftw_execute(plan_r2c_f);

    for (int i = 0; i < Nx * K; ++i) {
        double a = fft_f[i][0], b = fft_f[i][1];
        double c = fft_g_input[i][0], d = fft_g_input[i][1];

        fft_h[i][0] = a*c - b*d;
        fft_h[i][1] = a*d + b*c;
    }

    fftw_execute(plan_c2r);

    for (int i = 0; i < N; i++)
        h2[i] = (fft_out_real[i] / N);
}



void add_c(double *c1, double *c1_temp)
{
	for(int i=0;i<N;i++)
	{
		c1[i]-=c1_temp[i];
	}
}

void add_cv(double *c1, double *c1_temp)
{
	for(int i=0;i<Nx*Nz;i++)
	{
		c1[i]+=c1_temp[i];
	}
}



void getc1_fmt()
{
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			{
				dphidn0[i*Nz+j] = -log(1-n3[i*Nz+j]);
				dphidn1[i*Nz+j] =  n2[i*Nz+j]/(1-n3[i*Nz+j]);
				dphidn2[i*Nz+j] = n1[i*Nz+j]/(1-n3[i*Nz+j]) + (3*n2[i*Nz+j]*n2[i*Nz+j] - 3*(n2vx[i*Nz+j]*n2vx[i*Nz+j] + n2vz[i*Nz+j]*n2vz[i*Nz+j]))/(24*PI*(1-n3[i*Nz+j])*(1-n3[i*Nz+j]));
				dphidn3[i*Nz+j] = n0[i*Nz+j]/(1-n3[i*Nz+j]) + (n1[i*Nz+j]*n2[i*Nz+j] - n1vx[i*Nz+j]*n2vx[i*Nz+j] - n1vz[i*Nz+j]*n2vz[i*Nz+j])/(1-n3[i*Nz+j])/(1-n3[i*Nz+j]) + (n2[i*Nz+j]*n2[i*Nz+j]*n2[i*Nz+j] - 3*n2[i*Nz+j]*(n2vx[i*Nz+j]*n2vx[i*Nz+j] + n2vz[i*Nz+j]*n2vz[i*Nz+j]))/12/PI/(1-n3[i*Nz+j])/(1-n3[i*Nz+j])/(1-n3[i*Nz+j]);
				dphidn1vx[i*Nz+j] = -n2vx[i*Nz+j]/(1-n3[i*Nz+j]);
				dphidn1vz[i*Nz+j] = -n2vz[i*Nz+j]/(1-n3[i*Nz+j]);
				dphidn2vx[i*Nz+j] = -n1vx[i*Nz+j]/(1-n3[i*Nz+j]) - 3*n2[i*Nz+j]*n2vx[i*Nz+j]/(12*PI*(1-n3[i*Nz+j])*(1-n3[i*Nz+j]));
				dphidn2vz[i*Nz+j] = -n1vz[i*Nz+j]/(1-n3[i*Nz+j]) - 3*n2[i*Nz+j]*n2vz[i*Nz+j]/(12*PI*(1-n3[i*Nz+j])*(1-n3[i*Nz+j]));
				
			}
	
	for(int i1=0;i1<Nx;i1++)
		for(int j1=0;j1<Nz;j1++)
			{
				c1[0][i1*Nz+j1]=0.;
				c1[1][i1*Nz+j1]=0.;
				c1[2][i1*Nz+j1]=0.;
				//c1_temp[i*Nz+j1]=0;
			}
	
	conv_FFT2D_2(dphidn0, omega0, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
	conv_FFT2D_2(dphidn1, omega1, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
	
	conv_FFT2D_2(dphidn2, omega2, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
	
	conv_FFT2D_2(dphidn3, omega3, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
	
	conv_FFT2D_2(dphidn1vx, omega1vx, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
	
	conv_FFT2D_2(dphidn1vz, omega1vz, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
	
	conv_FFT2D_2(dphidn2vx, omega2vx, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
	
	conv_FFT2D_2(dphidn2vz, omega2vz, c1_temp);
	add_c(c1[0],c1_temp);
	add_c(c1[1],c1_temp);
	add_c(c1[2],c1_temp);
	
				
}


void getVext(){

	for(int i=0;i<N;i++)
		Vext[i]=0.0;
	
	for(int i=0;i<Nx;i++)
		for(int j=0;j<NiR-1;j++)
			Vext[i*Nz+j] = 1000.;
	
	
	for(int i=0;i<Nx;i++)
	for(int j=NiR-1;j<Nz;j++)
			Vext[i*Nz+j] = eps*ew*(2./15*pow(dz*j,-9)-pow(dz*j,-3));
	
			
	for(int i=0;i<Nx;i++)
	{
		psi_bottom[i]=Vq;
		psi_top[i]=0.;
	}
	
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz;j++)
			Vext_Q[i*Nz+j] = 0.;//(dz*(Nz/2-j))/Lz*.2*(8.85e-12);
}


void getc1_LJ2()
{

	getc1_LJ2_i(0,0,1.);
	getc1_LJ2_i(0,1,eps_i/eps);
	getc1_LJ2_i(0,2,eps_i/eps);
	
	getc1_LJ2_i(1,0,eps_i/eps);
	getc1_LJ2_i(2,0,eps_i/eps);

	
}

void getc1_LJ2_i(int j1 ,int k1 ,double factor)
{
	for(int i=0;i<Nx;i++)
		for(int j=0;j<Nz-NiLJ+1;j++)
			rhocopy[i*Nz+j]=rho[k1][i*Nz+j]*factor;
			
	for(int i=0;i<Nx;i++)
		for(int j=Nz-NiLJ+1;j<Nz;j++)
			rhocopy[i*Nz+j]=0.0;
	
	conv_FFT2D_2(rhocopy,UfilterFFT, c1_temp);
	add_c(c1[j1],c1_temp);
}


void getUfilterFFT()
{
	FILE *F=fopen("../data/U2k.dat","r");
	
	double k,Fk,a,b,Fb;
	
	double Utempf[1000];
	int i=0;
	while(fscanf(F,"%lf %lf %lf %lf %lf",&a,&k,&b,&Fk,&Fb)==5)
	{
		Utempf[i++]=Fk;
	}
	fclose(F);
	Utempf[0]=Fb;
	
	
	double kx=0.0,kz=0.0;
  double dkx=2.*PI/Lx,dkz=2.*PI/Lz;
  for(int i=0;i<Nx;i++)
  {
  	if(i<=Nx/2)
  		kx=dkx*i;
  	else
  		kx=-dkx*(Nx-i);
  	for(int j=0;j<K;j++)
  	{
  		kz=dkz*j;
			k=sqrt(kx*kx+kz*kz);
			if(i==0 && j==0)
			{
				UfilterFFT[i*K+j][0]=-eps*14.726050906548;
				UfilterFFT[i*K+j][1]=0.0;
			}
			else
			{
				int key=MIN((int)(k/PI*100),1000-1);
				
				UfilterFFT[i*K+j][0]=eps*Utempf[key]/2;
				//printf("%lf %lf\n",k,UfilterFFT[i*K+j][0]);
				UfilterFFT[i*K+j][1]=0.0;
  		}
  	}
  }
	
}


void rhocpy()
{
	for(int i=0;i<N;i++)
		rhocopy[i]=(rho[0][i]+rho[1][i]+rho[2][i]);
	
	for(int i=0;i<Nx;i++)
		for(int j=Nz-NiR+1;j<Nz;j++)
			 rhocopy[i*Nz+j]= 0.0 ;
}

void filterrho()
{
	for(int i=0;i<Nx;i++)
		for(int j=Nz-3*NiLJ;j<Nz;j++)
			{
				rho[0][i*Nz+j]=rhob_solvent_g;
				rho[1][i*Nz+j]=rhob_solute_g;
				rho[2][i*Nz+j]=rhob_solute_g;
			}
}


void getpsi()
{
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nz;j++)
		phi[j*Nx+i]=-(rho[1][i*Nz+j]-rho[2][i*Nz+j]);
	}
	
	poisson_2D_complex(dx, Nx, Nz,psi_bottom,psi_top,phi,psi);
}

void iterate(){	
	rhocpy();
	getn();
	getc1_fmt();
	getc1_LJ2();
	
	
	getpsi();
	
	double mu[3]={mu_solvent,mu_solute,mu_solute};
	
	int iend=Nz-3*NiLJ;
	
	
	for(int i=0;i<Nx;i++)
		for(int j=NiR;j<iend;j++)
			{
				rhonew[0][i*Nz+j]=exp(-Vext[i*Nz+j]+c1[0][i*Nz+j]+mu[0]);
				rhonew[1][i*Nz+j]=exp(-Vext[i*Nz+j]+c1[1][i*Nz+j]-lambdaB*psi[j*Nx+i]+mu[1]-Vext_Q[i*Nz+j]);
				rhonew[2][i*Nz+j]=exp(-Vext[i*Nz+j]+c1[2][i*Nz+j]+lambdaB*psi[j*Nx+i]+mu[2]+Vext_Q[i*Nz+j]);
			}
			
			
	for(int i=0;i<Nx;i++)
		for(int j=NiR;j<iend;j++)
			{
				rho[0][i*Nz+j]=alpha*rhonew[0][i*Nz+j] + (1-alpha)*rho[0][i*Nz+j];
				rho[1][i*Nz+j]=alpha*rhonew[1][i*Nz+j] + (1-alpha)*rho[1][i*Nz+j];
				rho[2][i*Nz+j]=alpha*rhonew[2][i*Nz+j] + (1-alpha)*rho[2][i*Nz+j];
			}

	double Nrhotemp[3]={0.0,0.0,0.0};
	for(int i=0;i<Nx;i++)
	for(int j=NiR;j<Nz;j++)
		{
			Nrhotemp[0]+=(rho[0][i*Nz+j]-rhob_solvent_g);
			Nrhotemp[1]+=(rho[1][i*Nz+j]-rhob_solute_g);
			Nrhotemp[2]+=(rho[2][i*Nz+j]-rhob_solute_g);
		}
		
	for(int j=0;j<3;j++)
		Nrhotemp[j]=(((Nrhotemp[j]*dx*dz)/Lx));
	
	for(int i=0;i<Nx;i++)
	for(int j=NiR;j<Nz;j++)
		{
			rho[0][i*Nz+j]=(rho[0][i*Nz+j]-rhob_solvent_g)*(Nrho[0]/Nrhotemp[0])+rhob_solvent_g;
			rho[1][i*Nz+j]=(rho[1][i*Nz+j]-rhob_solute_g)*(Nrho[1]/Nrhotemp[1])+rhob_solute_g;
			rho[2][i*Nz+j]=(rho[2][i*Nz+j]-rhob_solute_g)*(Nrho[2]/Nrhotemp[2])+rhob_solute_g;
		}


}


void write_rho(double elapsed,int count_iter)
{
	char fname[500];
	sprintf(fname,"../data/newrho2DNeps%feps_i%few%fdx%frhob_solvent_l%frhob_solute_l%flambdaB%fVq%f.dat",eps,eps_i,ew,dx,rhob_solvent_l,rhob_solute_l,lambdaB,Vq);
	FILE *F=fopen(fname,"w");
	for(int i=0;i<Nx;i++)
	{
		for(int j=0;j<Nz;j++)
			fprintf(F,"%d %d %f %f %f %f\n",i,j,rho[0][i*Nz+j],rho[1][i*Nz+j],rho[2][i*Nz+j],psi[j*Nx+i]);
		fprintf(F,"\n");
	}
	fprintf(F,"------------------\ncycles: %d x %d ;time: %f s\n",Nbatch,count_iter,elapsed);
	fclose(F);
	
	sprintf(fname,"../data/newrho1DNeps%feps_i%few%fdx%frhob_solvent_l%frhob_solute_l%flambdaB%fVq%f.dat",eps,eps_i,ew,dx,rhob_solvent_l,rhob_solute_l,lambdaB,Vq);
	F=fopen(fname,"w");
	
	for(int j=0;j<Nz;j++)
	fprintf(F,"%f %f %f %f %f %f\n",j*dz,rho[0][(Nx/2)*Nz+j],rho[1][(Nx/2)*Nz+j],rho[2][(Nx/2)*Nz+j],c1[0][(Nx/2)*Nz+j],psi[j*Nx+(Nx/2)]);
	fprintf(F,"------------------\ncycles: %d x %d ;time: %f s\n",Nbatch,count_iter,elapsed);
	fclose(F);
}

void initialize_dataframes(int Narr,double *arr)
{
	for(int i=0;i<Narr;i++)
	arr[i]=0.0;
}

void initialize_all_dataframes()
{
	for(int i=0;i<3;i++)
	{
		initialize_dataframes(N,rho[i]);
		initialize_dataframes(N,rhonew[i]);
		initialize_dataframes(N,c1[i]);
	}
	
	initialize_dataframes(N,rhocopy);
	initialize_dataframes(N,Vext);
	initialize_dataframes(N,Vext_Q);
	initialize_dataframes(N,phi);
	initialize_dataframes(N,psi);
	initialize_dataframes(Nx,psi_bottom);
	initialize_dataframes(Nx,psi_top);
	
	initialize_dataframes(N,n0);
	initialize_dataframes(N,n1);
	initialize_dataframes(N,n2);
	initialize_dataframes(N,n3);
	initialize_dataframes(N,n1vx);
	initialize_dataframes(N,n1vz);
	initialize_dataframes(N,n2vx);
	initialize_dataframes(N,n2vz);
	
	
	initialize_dataframes(N,dphidn0);
	initialize_dataframes(N,dphidn1);
	initialize_dataframes(N,dphidn2);
	initialize_dataframes(N,dphidn3);
	initialize_dataframes(N,dphidn1vx);
	initialize_dataframes(N,dphidn1vz);
	initialize_dataframes(N,dphidn2vx);
	initialize_dataframes(N,dphidn2vz);
	
	initialize_dataframes(N,c1_temp);
}

void init_fftw()
{
    fft_in  = (double*) fftw_malloc(sizeof(double) * N);
    fft_in2 = (double*) fftw_malloc(sizeof(double) * N);
    fft_out_real = (double*) fftw_malloc(sizeof(double) * N);

    fft_f = fftw_alloc_complex((size_t)Nx * K);
    fft_g = fftw_alloc_complex((size_t)Nx * K);
    fft_h = fftw_alloc_complex((size_t)Nx * K);

    plan_r2c_f = fftw_plan_dft_r2c_2d(Nx, Nz, fft_in, fft_f, FFTW_MEASURE);
    plan_r2c_g = fftw_plan_dft_r2c_2d(Nx, Nz, fft_in2, fft_g, FFTW_MEASURE);
    plan_c2r   = fftw_plan_dft_c2r_2d(Nx, Nz, fft_h, fft_out_real, FFTW_MEASURE);
}

void cleanup_fftw()
{
    fftw_destroy_plan(plan_r2c_f);
    fftw_destroy_plan(plan_r2c_g);
    fftw_destroy_plan(plan_c2r);

    fftw_free(fft_in);
    fftw_free(fft_in2);
    fftw_free(fft_out_real);

    fftw_free(fft_f);
    fftw_free(fft_g);
    fftw_free(fft_h);

    fftw_cleanup();
}

int main(int argc,char *argv[])
{
	
	read_params();
	initialize_all_dataframes();
	init_fftw();
	
	clock_t start = clock();clock_t end;
	double elapsed;
	
	get_mu();
	getVext();
	rhoinit();
	getUfilterFFT();
	getomega3();
	getomega2();
	
	getomega1();
	getomega0();
	getomega2v();
	getomega1v();
  
	
	for(int i=1;i<=INT_MAX;++i)
	{	
		for(int j=1;j<=Nbatch;j++)
			iterate();
		end = clock();
		elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
		write_rho(elapsed,i);
	}
	cleanup_fftw();

	return 0;
}
