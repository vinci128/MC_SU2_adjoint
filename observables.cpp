#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "ranlxd.h"
#include "complex.h"
#include "global.h"
#include "geometry.h"
#include "update.h"
#include "su2.h"
#include "observables.h"

double plaquette(int s, int mu,int nu){

double plaq=0;

int sp[Nd];


// first neighbor of the site s in direction nu

	sp[nu] = neig[s][nu];
	sp[mu] = neig[s][mu];

su2_x Unu,Umu,Unu_adj,Umu_adj;
su2_x temp1,temp2,prod;

Unu = U[s][nu];
Umu = U[sp[nu]][mu];


_su2_adj(Unu_adj,Unu);
_su2_adj(Umu_adj,Umu);


_su2_times_su2(temp1,Umu_adj,Unu_adj);
_su2_times_su2(temp2,U[sp[mu]][nu],temp1);
_su2_times_su2(prod,U[s][mu],temp2);
// TrU = x_0 *Tr id = 2 x_0
plaq = 2*prod.x[0];

return plaq;
}

double plaq_smear(int s, int mu,int nu){

double plaq=0;

int sp[Nd];


// first neighbor of the site s in direction nu

	sp[nu] = neig[s][nu];
	sp[mu] = neig[s][mu];

su2_x Unu,Umu,Unu_adj,Umu_adj;
su2_x temp1,temp2,prod;

Unu = U_smear[s][nu];
Umu = U_smear[sp[nu]][mu];


_su2_adj(Unu_adj,Unu);
_su2_adj(Umu_adj,Umu);


_su2_times_su2(temp1,Umu_adj,Unu_adj);
_su2_times_su2(temp2,U_smear[sp[mu]][nu],temp1);
_su2_times_su2(prod,U_smear[s][mu],temp2);
// TrU = x_0 *Tr id = 2 x_0
plaq = 2*prod.x[0];

return plaq;
}

double avr_plaquette(){

double plaq_average=0;
double plaq=0;

for (int s = 0; s < V; s++) {
for (int nu = 0; nu < Nd; nu++) {
for( int mu =0 ; mu < nu; mu++){

plaq = plaquette(s,mu,nu);
plaq_average += plaq;

}
}
}



plaq_average = plaq_average/(12*V);

return plaq_average;

}

double avr_plaquette_smear(){

double plaq_average=0;
double plaq=0;

for (int s = 0; s < V; s++) {
for (int nu = 0; nu < Nd; nu++) {
for( int mu =0 ; mu < nu; mu++){

plaq = plaq_smear(s,mu,nu);
plaq_average += plaq;

}
}
}



plaq_average = plaq_average/(12*V);

return plaq_average;

}

double phi_sq(){

double phi_average=0;

for (int s = 0; s < V; s++) {
   for(int a = 0; a < Ng; a++){

      phi_average+=phi[s][a]*phi[s][a];

      }
}

phi_average = phi_average/(V);

return phi_average;

}

double phi_sq_smear(){

double phi_average=0;

for (int s = 0; s < V; s++) {
   for(int a = 0; a < Ng; a++){

      phi_average+=phi_smear[s][a]*phi_smear[s][a];

      }
}

phi_average = phi_average/(V);

return phi_average;

}

void SC1_t(double O[]){

	double Phi2;

			//printf(" Length at point %i: %1.8e\n", ix, Phi_sq );
			int t, x, y, z; for (t=0; t<Nt; t++) { O[t] = 0; for (z=0; z<Nz; z++) for (y=0; y<Ny; y++)  for (x=0; x<Nx; x++) {

				  int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;
				  Phi2=0;
				 for (int a = 0; a < Ng; ++a)
				 {
				 Phi2 += phi[s][a]*phi[s][a];
				 }
				O[t] += Phi2;
		}
		O[t] = O[t]/ sqrt(Nx*Ny*Nz);
	}

//printf("%f %f %f %f %f %f %f %f\n",OSC1[0], OSC1[1],OSC1[2],OSC1[3],OSC1[4],OSC1[5],OSC1[6],OSC1[7] );
//	return Phi_cond/(double)GLB_VOLUME;

}

void SC2_t(double O[]){

	double W[4][3][3];
	su2 U_x[4];
	su2 U_x_adj[4];

	su2 S[3];
	su2 S1,S2,S3;


	_su2_t1(S[0]);
	_su2_t2(S[1]);
	_su2_t3(S[2]);

	//double Phix[3], Phix_up[3];


	//double Hop = 0;
	int t, x, y, z; for (t=0; t<Nt; t++) { O[t] = 0; for (z=0; z<Nz; z++) for (y=0; y<Ny; y++)  for (x=0; x<Nx; x++) {

	  int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;

 int sp[Nd];

for(int mu= 0; mu < Nd; mu++){
sp[mu] = neig[s][mu];
}

		for (int mu = 0; mu < 3; mu++) {

			//Phix_up = *_FIELD_AT(u_adjoint, iup(ix,mu)); //Phi(x + mu)
			_su2_x_represent(U_x[mu],U[s][mu]); // U_mu(x)
			_suNg_dagger(U_x_adj[mu], U_x[mu]); //U_mu(x)^dag

			for (int a = 0; a < 3; a++) {
			 for (int b = 0; b < 3; b++) {

				 _suNg_times_suNg(S1,S[b],U_x_adj[mu]);
				 _suNg_times_suNg(S2,U_x[mu],S1);
				 _suNg_times_suNg(S3,S[a],S2);
				 _suNg_trace_re(W[mu][a][b],S3);

				 O[t] -= 2*phi[s][a]*W[mu][a][b]*phi[sp[mu]][b];

				}
			}
		}


	}
	O[t] = O[t]/ sqrt(3*Nx*Ny*Nz);
}

}

void SC3_t(double O[]){

	int t, x, y, z; for (t=0; t<Nt; t++) { O[t] = 0; for (z=0; z<Nz; z++) for (y=0; y<Ny; y++)  for (x=0; x<Nx; x++) {

		int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;


		O[t] += plaq_smear(s,0,1);
		O[t] += plaq_smear(s,1,2);
		O[t] += plaq_smear(s,2,0);

		}
	O[t] = O[t]/ sqrt(3*Nx*Ny*Nz);
	}

}

double phi_plaq_smear(int s,int mu,int nu)
{
  int sp[Nd];
  double p;
  su2 v1,v2,v3,v4,w1,w2,w3,w4;

for (int m = 0; m < Nd; m++)
{
	sp[m] = neig[s][m];
}

su2 Phim;
//printf("%f\n", Phix.c[2]);

(Phim).c[0].im = 0.0;
(Phim).c[0].re = +5.000000000000000e-01*phi_smear[s][2];
(Phim).c[1].im = -5.000000000000000e-01*phi_smear[s][1];
(Phim).c[1].re = +5.000000000000000e-01*phi_smear[s][0];
(Phim).c[2].im = 5.000000000000000e-01*phi_smear[s][1];
(Phim).c[2].re = +5.000000000000000e-01*phi_smear[s][0];
(Phim).c[3].im = 0.0;
(Phim).c[3].re = -5.000000000000000e-01*phi_smear[s][2];

_su2_x_represent(v1,U_smear[s][mu]);
_su2_x_represent(v2,U_smear[sp[mu]][nu]);
_su2_x_represent(v3,U_smear[sp[nu]][mu]);
_su2_x_represent(v4,U_smear[s][nu]);

  _suNg_times_suNg(w1,v1,v2);
  _suNg_times_suNg(w2,v4,v3);
  _suNg_times_suNg_dagger(w3,w1,w2);
	_suNg_times_suNg(w4,Phim,w3);

_suNg_trace_im(p,w4);

//printf("Phi_plaq_im = %f\n", p);
  return p;
}

void B_t(double O[][3]){

	//suNg_algebra_vector Phic;
	//double corr = 0;
	//double SC1 = 0;
	//double O[8] ={0};
	double phi2;
	double phip[3];

			//printf(" Length at point %i: %1.8e\n", ix, Phi_sq );
	int t, x, y, z;
	for (t=0; t<Nt; t++) {

			O[t][0] = 0;
			O[t][1] = 0;
			O[t][2] = 0;

			for (z=0; z<Nz; z++) {
			for (y=0; y<Ny; y++) {
			for (x=0; x<Nx; x++) {

				int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;
				phi2=0;
				for (int a = 0; a < Ng; ++a)
					{
						phi2 += phi_smear[s][a]*phi_smear[s][a];
					}


				phip[0] = phi_plaq_smear(s,1,2);
				phip[1] = phi_plaq_smear(s,2,0);
				phip[2] = phi_plaq_smear(s,0,1);

				for (int i = 0; i < 3; i++) {
					O[t][i] += phip[i]/sqrt(Nx*Ny*Nz*phi2);
				}
		}}}
	}

}

void B_z_t(double O[][3]){

	//suNg_algebra_vector Phic;
	//double corr = 0;
	//double SC1 = 0;
	//double O[8] ={0};
	double phi2;
	double phip[3];
	vec p_z;
	double omega;
p_z.v[0]=0;
p_z.v[1]=0;
p_z.v[2]=2*M_PI/Nz;
			//printf(" Length at point %i: %1.8e\n", ix, Phi_sq );
vec r ;
	int t, x, y, z;
	for (t=0; t<Nt; t++) {

O[t][0]=0;
O[t][1]=0;
O[t][2]=0;


			for (z=0; z<Nz; z++) {
			for (y=0; y<Ny; y++) {
			for (x=0; x<Nx; x++) {
r.v[0] =x;
r.v[1] =y;
r.v[2] =z;

_sc_prod(omega,r,p_z);

				int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;
				phi2=0;
for (int a = 0; a < Ng; ++a)
{
phi2 += phi_smear[s][a]*phi_smear[s][a];
}


				phip[0] = phi_plaq_smear(s,1,2);
				phip[1] = phi_plaq_smear(s,2,0);
				phip[2] = phi_plaq_smear(s,0,1);

				for (int i = 0; i < 3; i++) {
					O[t][i] += phip[i]/sqrt(Nx*Ny*Nz*phi2)*cos(omega);


				}
		}}}
	}

}

void B2_z_t(double O[][3]){

	//suNg_algebra_vector Phic;
	//double corr = 0;
	//double SC1 = 0;
	//double O[8] ={0};
	double phi2;
	double phip[3];
	vec p_z;
	double omega;
p_z.v[0]=0;
p_z.v[1]=0;
p_z.v[2]=2*M_PI/Nz;
			//printf(" Length at point %i: %1.8e\n", ix, Phi_sq );
vec r ;
	int t, x, y, z;
	for (t=0; t<Nt; t++) {

O[t][0]=0;
O[t][1]=0;
O[t][2]=0;


			for (z=0; z<Nz; z++) {
			for (y=0; y<Ny; y++) {
			for (x=0; x<Nx; x++) {
r.v[0] =x;
r.v[1] =y;
r.v[2] =z;

_sc_prod(omega,r,p_z);

				int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;
				phi2=0;
for (int a = 0; a < Ng; ++a)
{
phi2 += phi_smear[s][a]*phi_smear[s][a];
}


				phip[0] = phi_plaq_smear(s,1,2)*phi_plaq_smear(s,1,2);
				phip[1] = phi_plaq_smear(s,2,0)*phi_plaq_smear(s,2,0);
				phip[2] = phi_plaq_smear(s,0,1)*phi_plaq_smear(s,0,1);

				for (int i = 0; i < 3; i++) {
					O[t][i] += phip[i]/sqrt(Nx*Ny*Nz*phi2)*cos(omega);


				}
		}}}
	}

}

void Bphi_z_t(double O[][3]){

	//suNg_algebra_vector Phic;
	//double corr = 0;
	//double SC1 = 0;
	//double O[8] ={0};
	double phi2;
	double phip[3];
	vec p_z;
	double omega;
p_z.v[0]=0;
p_z.v[1]=0;
p_z.v[2]=2*M_PI/Nz;
			//printf(" Length at point %i: %1.8e\n", ix, Phi_sq );
vec r ;
	int t, x, y, z;
	for (t=0; t<Nt; t++) {

O[t][0]=0;
O[t][1]=0;
O[t][2]=0;


			for (z=0; z<Nz; z++) {
			for (y=0; y<Ny; y++) {
			for (x=0; x<Nx; x++) {
r.v[0] =x;
r.v[1] =y;
r.v[2] =z;

_sc_prod(omega,r,p_z);

				int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;
				phi2=0;
for (int a = 0; a < Ng; ++a)
{
phi2 += phi_smear[s][a]*phi_smear[s][a];
}


				phip[0] = phi2*phi_plaq_smear(s,1,2);
				phip[1] = phi2*phi_plaq_smear(s,2,0);
				phip[2] = phi2*phi_plaq_smear(s,0,1);

				for (int i = 0; i < 3; i++) {
					O[t][i] += phip[i]/sqrt(Nx*Ny*Nz*phi2)*cos(omega);


				}
		}}}
	}

}

void B_2z_t(double O[][3]){

	//suNg_algebra_vector Phic;
	//double corr = 0;
	//double SC1 = 0;
	//double O[8] ={0};
	double phi2;
	double phip[3];
	vec p_z;
	double omega;
p_z.v[0]=0;
p_z.v[1]=0;
p_z.v[2]=4*M_PI/Nz;
			//printf(" Length at point %i: %1.8e\n", ix, Phi_sq );
vec r ;
	int t, x, y, z;
	for (t=0; t<Nt; t++) {

O[t][0]=0;
O[t][1]=0;
O[t][2]=0;


			for (z=0; z<Nz; z++) {
			for (y=0; y<Ny; y++) {
			for (x=0; x<Nx; x++) {
r.v[0] =x;
r.v[1] =y;
r.v[2] =z;

_sc_prod(omega,r,p_z);

				int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;
				phi2=0;
for (int a = 0; a < Ng; ++a)
{
phi2 += phi_smear[s][a]*phi_smear[s][a];
}


				phip[0] = phi_plaq_smear(s,1,2);
				phip[1] = phi_plaq_smear(s,2,0);
				phip[2] = phi_plaq_smear(s,0,1);

				for (int i = 0; i < 3; i++) {
					O[t][i] += phip[i]/sqrt(Nx*Ny*Nz*phi2)*cos(omega);


				}
		}}}
	}

}


void W_t(double O[][6][6][6]){

int sp[Nd], sm[Nd], spp[Nd][Nd],spm[Nd][Nd],sppm[Nd][Nd][Nd];

int t, x, y, z;
for (t=0; t<Nt; t++) {

		for (z=0; z<Nz; z++) {
		for (y=0; y<Ny; y++) {
		for (x=0; x<Nx; x++) {

int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;

// Initialize neighbours and observables arrays

	for (int mu = 0; mu < Nd; mu++) {
sp[mu] = neig[s][mu];
sm[mu] = neig[s][Nd + mu];
		for (int nu = 0; nu < Nd; nu++) {
			spp[mu][nu] = neig[sp[mu]][nu];
			spm[mu][nu] = neig[sp[mu]][Nd+nu];
				for (int rho = 0; rho < Nd; rho++) {
					sppm[mu][nu][rho] = neig[spp[mu][nu]][Nd+rho];

				}
			}
}
					for (int mu = 0; mu < 3; mu++) {
						for (int nu = 0; nu < 3; nu++) {
								for (int rho = 0; rho < 3; rho++) {

									O[t][mu][nu][rho] = 0;
									O[t][mu+3][nu+3][rho+3] = 0;

				su2_x v1,v2,v3,v4,v5,v6,v7,v8;
				su2_x w1,w2,w3,w4,w5,w6,w7,w8;
				su2_x S1,S2,S3,S4,S5,S6,S7;
				//double W,W_min;

			v1=U_smear[s][mu];
			v2=U_smear[sp[mu]][mu];
			v3=U_smear[spp[mu][mu]][nu];
			//v4=*pu_gauge_smear(im,mu);
			_su2_adj(v4,U_smear[spp[mu][nu]][mu]);
			v5=U_smear[spp[mu][nu]][rho];
			_su2_adj(v6,U_smear[spp[nu][rho]][mu]);
			_su2_adj(v7,U_smear[sp[nu]][rho]);
			_su2_adj(v8,U_smear[s][nu]);



			_su2_times_su2(S1,v7,v8);
			_su2_times_su2(S2,v6,S1);
			_su2_times_su2(S3,v5,S2);
			_su2_times_su2(S4,v4,S3);
			_su2_times_su2(S5,v3,S4);
			_su2_times_su2(S6,v2,S5);
			_su2_times_su2(S7,v1,S6);
			O[t][mu][nu][rho] += S7.x[0];

			_su2_adj(w1,U_smear[sm[mu]][mu]);
			_su2_adj(w2,U_smear[s][mu]);
			_su2_adj(w3,U_smear[sppm[mu][mu][nu]][nu]);
			w4=U_smear[sp[nu]][mu];
			_su2_adj(w5,U_smear[sppm[mu][nu][rho]][rho]);
			w6=U_smear[sppm[nu][rho][mu]][mu];
			w7=U_smear[spm[nu][rho]][rho];
			w8=U_smear[sm[nu]][nu];

			_su2_times_su2(S1,w7,w8);
			_su2_times_su2(S2,w6,S1);
			_su2_times_su2(S3,w5,S2);
			_su2_times_su2(S4,w4,S3);
			_su2_times_su2(S5,w3,S4);
			_su2_times_su2(S6,w2,S5);
			_su2_times_su2(S7,w1,S6);
			O[t][mu+3][nu+3][rho+3] += S7.x[0];


			}}}


					}
				}
			}

		}


}

void Bm_t(double O[][3][3][3] ){

double W[Nt][6][6][6];

W_t(W);

for (int t=0; t<Nt; t++) {

	for (int mu = 0; mu < 3; mu++) {
		for (int nu = 0; nu < 3; nu++) {
				for (int rho = 0; rho < 3; rho++) {

O[t][mu][nu][rho] = W[t][mu][nu][rho]+ W[t][mu][nu][rho+3] +W[t][mu][nu+3][rho] +W[t][mu][nu+3][rho+3]
 									- W[t][mu+3][nu][rho] -W[t][mu+3][nu][rho+3] -W[t][mu+3][nu+3][rho] - W[t][mu+3][nu+3][rho+3];


}
}
}

}

}

void Cm_t(double O[][3][3][3]){

	double W[Nt][6][6][6];

	W_t(W);

for (int t=0; t<Nt; t++) {

	for (int mu = 0; mu < 3; mu++) {
		for (int nu = 0; nu < 3; nu++) {
				for (int rho = 0; rho < 3; rho++) {

O[t][mu][nu][rho] = W[t][mu][nu][rho]+ W[t][mu][nu][rho+3] -W[t][mu][nu+3][rho] -W[t][mu][nu+3][rho+3]
					 				+ W[t][mu+3][nu][rho] +W[t][mu+3][nu][rho+3] -W[t][mu+3][nu+3][rho] - W[t][mu+3][nu+3][rho+3];


}
}
}

}

}

void Dm_t(double O[][3][3][3]){

double W[Nt][6][6][6];
W_t(W);

for (int t=0; t<Nt; t++) {

	for (int mu = 0; mu < 3; mu++) {
		for (int nu = 0; nu < 3; nu++) {
				for (int rho = 0; rho < 3; rho++) {

O[t][mu][nu][rho] = W[t][mu][nu][rho]- W[t][mu][nu][rho+3] +W[t][mu][nu+3][rho] -W[t][mu][nu+3][rho+3]
									+ W[t][mu+3][nu][rho] -W[t][mu+3][nu][rho+3] +W[t][mu+3][nu+3][rho] - W[t][mu+3][nu+3][rho+3];


}
}
}

}

}

void T1m_t(double O[][3]){

double B[Nt][3][3][3];

Bm_t(B);

for (int t=0; t<Nt; t++) {

	O[t][0] = B[t][0][1][2]+B[t][0][2][1];
	O[t][1] = B[t][1][2][0]+B[t][1][0][2];
	O[t][2] = B[t][2][0][1]+B[t][2][1][0];
}

}

void T2m_t(double O[][3]){

	double C[Nt][3][3][3];

	Cm_t(C);

	for (int t=0; t<Nt; t++) {

		O[t][0] = C[t][0][1][2]+C[t][2][1][0];
		O[t][1] = C[t][1][2][0]+C[t][0][2][1];
		O[t][2] = C[t][2][0][1]+C[t][1][0][2];

	}

}

void T3m_t(double O[][3]){

	double D[Nt][3][3][3];

	Dm_t(D);

	for (int t=0; t<Nt; t++) {

		O[t][0] = D[t][0][1][2]+D[t][1][0][2];
		O[t][1] = D[t][1][2][0]+D[t][2][1][0];
		O[t][2] = D[t][2][0][1]+D[t][0][2][1];

	}

}
