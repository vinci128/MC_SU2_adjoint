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
  ;
  double p;
  su2 v1,v2,v3,v4,w1,w2,w3,w4;

for (int mu = 0; mu < Nd; mu++)
{
	sp[mu] = neig[s][mu];
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
//printf("%f %f %f %f %f %f %f %f\n",OSC1[0], OSC1[1],OSC1[2],OSC1[3],OSC1[4],OSC1[5],OSC1[6],OSC1[7] );
//	return Phi_cond/(double)GLB_VOLUME;

}
