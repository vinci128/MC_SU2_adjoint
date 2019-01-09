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

#define M_PI 3.14159265358979323846

double accepted =0;
double op=0;
double eta;
double eps = 0.2;

double lambda ;
double kappa ;
double beta ;


 void cold_start(){


// Initialization of the fields U_mu(x) = x_0 id + i x_i sigma_i = Id, phi(x) = phi^a(x) T^a = 0

  for(int s=0; s < V; s++){

      for(int mu =0; mu < Nd; mu++){
    U[s][mu].x[0]=1;
    U[s][mu].x[1]=0;
    U[s][mu].x[2]=0;
    U[s][mu].x[3]=0;
 //     std::cout << "Field at point: " << U[s][mu].x[0] << "\n";
     }

          for(int a =0; a < Ng; a++){
    phi[s][a]=0;
     }

  }

}

 void hot_start(){
// Initialization of the fields U_mu(x) = x_0 id + i x_i sigma_i , phi(x) = phi^a(x) T^a

  for(int s=0; s < V; s++){

    double r[4];
ranlxd(r,4);

 double p[3];
ranlxd(p,3);

double N =1/sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2] +r[3]*r[3]);
r[0]= r[0]*N;
r[1]= r[1]*N;
r[2]= r[2]*N;
r[3]= r[3]*N;

      for(int mu =0; mu < Nd; mu++){
    U[s][mu].x[0]=r[0];
    U[s][mu].x[1]=r[1];
    U[s][mu].x[2]=r[2];
    U[s][mu].x[3]=r[3];
 //     std::cout << "Field at point: " << U[s][mu].x[0] << "\n";
     }

          for(int a =0; a < Ng; a++){
    phi[s][a]=2*p[a]-1;
     }

  }

}


void gauge_update(int s, int mu,double eps){

double t[5];
ranlxd(t,5);
double modr;
double r0 = t[0]- 0.5;
double r1 = t[1] - 0.5;
double r2 = t[2] - 0.5;
double r3 = t[3] - 0.5;
//printf("r0: %f, r1: %f, r2: %f, r3: %f \n",r0, r1, r2, r3 );
//printf("%f\n",  r[  Nd*mu+ Nd*Nd*x+ Nd*Nx*Ny*y+ Nd*Nx*Ny*Nz*z+ Nd*Nx*Ny*Nz*Nt*t ]);
modr= sqrt(r1*r1 + r2*r2+ r3*r3);


X[s][mu].x[0] = r0/(fabs(r0))*sqrt(1.- eps*eps);
//printf("x0: %f\n", r0/fabs(r0));
X[s][mu].x[1] = eps*r1/modr;
X[s][mu].x[2] = eps*r2/modr;
X[s][mu].x[3] = eps*r3/modr;

if(t[4]<0.5){
for (int k = 1; k < 4; ++k)
{
  X[s][mu].x[k] = - X[s][mu].x[k];
}
}


_su2_times_su2(U_new[s][mu],X[s][mu],U[s][mu]);

}

su2_x staple(int s, int mu){

  int sp[Nd];
int sm[Nd];
int spm[Nd][Nd];

for(int mu= 0; mu < Nd; mu++){
sp[mu] = neig[s][mu];
sm[mu] = neig[s][Nd+mu];
for(int nu= 0; nu < Nd; nu++){
spm[mu][nu] = neig[neig[s][mu]][Nd+nu];
}
}

su2_x St;

_su2_null(St);

// we calculate here the sum of the staples

for(int nu =0 ; nu < Nd; nu++){
if(nu == mu) continue;
//logf << "a " << mu << std::endl;
  su2_x Unu,Umu,Unu_adj,Umu_adj,Upmu;
  su2_x U_pmumnu, Upm_adj;
  su2_x U_mnu, Umnu_adj;
  su2_x temp1,temp2,Stup,Stdown;

// Calculation of the "up" staple

  Unu = U[s][nu]; //  U(x)_nu
//  printf(" Trace of the link at site %d in direction %d : %f\n",s,nu, 2*Unu.x[0]);

  Umu = U[sp[nu]][mu];  // U(x+ nu)_mu

    Upmu = U[sp[mu]][nu]; // U(x+mu)_nu

  _su2_adj(Unu_adj,Unu); // U(x)_nu^dag
  _su2_adj(Umu_adj,Umu); // U(x+nu)_mu^dag

  _su2_times_su2(temp1,Umu_adj,Unu_adj); // U(x+nu)_mu^dag U(x)_nu^dag
  _su2_times_su2(Stup,Upmu,temp1); // U(x+mu)_nu U(x+nu)_mu^dag U(x)_nu^dag

// Calculation of the "down" staple

  U_pmumnu = U[spm[mu][nu]][nu]; // U(x+ mu - nu)_nu

_su2_adj(Upm_adj,U_pmumnu);  // U(x+ mu - nu)_nu^dag

 U_mnu = U[sm[nu]][mu];   //  U(x -nu)_mu
_su2_adj(Umnu_adj, U_mnu);   //  U(x -nu)_mu^dag

_su2_times_su2(temp2,Umnu_adj,U[sm[nu]][nu]);  //  U(x-nu)_mu^\dag  U(x-nu)_nu
_su2_times_su2(Stdown,Upm_adj,temp2); // U(x+ mu -nu)_nu^\dag  U(x-nu)_mu^\dag  U(x-nu)_nu

_su2_add(St,Stup,Stdown);


}

return St;

}

su2_x staple_in(su2_x **in,int s, int mu){

  int sp[Nd];
int sm[Nd];
int spm[Nd][Nd];

for(int mu= 0; mu < Nd; mu++){
sp[mu] = neig[s][mu];
sm[mu] = neig[s][Nd+mu];
for(int nu= 0; nu < Nd; nu++){
spm[mu][nu] = neig[neig[s][mu]][Nd+nu];
}
}

su2_x St;

_su2_null(St);

// we calculate here the sum of the staples

for(int nu =0 ; nu < 3; nu++){
if(nu == mu) continue;
//logf << "a " << mu << std::endl;
  su2_x Unu,Umu,Unu_adj,Umu_adj,Upmu;
  su2_x U_pmumnu, Upm_adj;
  su2_x U_mnu, Umnu_adj;
  su2_x temp1,temp2,Stup,Stdown;

// Calculation of the "up" staple

  Unu = in[s][nu]; //  U(x)_nu
//  printf(" Trace of the link at site %d in direction %d : %f\n",s,nu, 2*Unu.x[0]);

  Umu = in[sp[nu]][mu];  // U(x+ nu)_mu

    Upmu = in[sp[mu]][nu]; // U(x+mu)_nu

  _su2_adj(Unu_adj,Unu); // U(x)_nu^dag
  _su2_adj(Umu_adj,Umu); // U(x+nu)_mu^dag

  _su2_times_su2(temp1,Umu_adj,Unu_adj); // U(x+nu)_mu^dag U(x)_nu^dag
  _su2_times_su2(Stup,Upmu,temp1); // U(x+mu)_nu U(x+nu)_mu^dag U(x)_nu^dag

// Calculation of the "down" staple

  U_pmumnu = in[spm[mu][nu]][nu]; // U(x+ mu - nu)_nu

_su2_adj(Upm_adj,U_pmumnu);  // U(x+ mu - nu)_nu^dag

 U_mnu = in[sm[nu]][mu];   //  U(x -nu)_mu
_su2_adj(Umnu_adj, U_mnu);   //  U(x -nu)_mu^dag

_su2_times_su2(temp2,Umnu_adj,in[sm[nu]][nu]);  //  U(x-nu)_mu^\dag  U(x-nu)_nu
_su2_times_su2(Stdown,Upm_adj,temp2); // U(x+ mu -nu)_nu^\dag  U(x-nu)_mu^\dag  U(x-nu)_nu

_su2_add(St,Stup,Stdown);


}

return St;

}

double gauge_DeltaS(int s,int mu,su2_x St,double beta,double kappa){

  // We generate arrays of the neighbours to the site s

double delta =0;
int sp[Nd];

su2 Um, Um_new;

su2 T[3];

_su2_t1(T[0]);
_su2_t2(T[1]);
_su2_t3(T[2]);

  double TrPn, TrP;
su2_x Pnew,P;

//printf(" Trace of the staple : %f\n", 2*St.x[0]);

_su2_times_su2(Pnew,U_new[s][mu],St);
_su2_times_su2(P,U[s][mu],St);

TrPn = 2*Pnew.x[0];
//printf(" Trace of new plaquette : %f\n", TrPn);
TrP = 2*P.x[0];
//printf(" Trace of old plaquette : %f\n", TrP);

delta -= beta/2*(TrPn- TrP);

sp[mu] = neig[s][mu];


_su2_x_represent(Um,U[s][mu]);
_su2_x_represent(Um_new,U_new[s][mu]);

for(int a=0; a< Ng; a++ ){
for(int b=0; b< Ng; b++ ){

double W, W_new;

su2 S1,S2,S3;

_suNg_times_suNg_dagger(S1,T[b],Um);
_suNg_times_suNg(S2,Um,S1);
_suNg_times_suNg(S3,T[a],S2);
_suNg_trace_re(W,S3);

_suNg_times_suNg_dagger(S1,T[b],Um_new);
_suNg_times_suNg(S2,Um_new,S1);
_suNg_times_suNg(S3,T[a],S2);
_suNg_trace_re(W_new,S3);

delta -= 2*kappa*phi[s][a]*(W_new-W)*phi[sp[mu]][b];

}
}

return delta;

}




int gauge_MC_test(int s, int mu,double DeltaS){
double r[1];
ranlxd(r,1);
double w,p;
  if(DeltaS<0 || DeltaS ==0 ){
_su2_eq(U[s][mu],U_new[s][mu]);
//logf << "Deltas: " << DeltaS << " Local update accepted at point: " << x << y << z << t << "\n";
//accepted ++;
//op ++;
return 1;
} if(DeltaS>0){
w=exp(-DeltaS);
p=r[0];
if(p<w){
_su2_eq(U[s][mu],U_new[s][mu]);
//accepted ++;
//op ++;
return 1;
//logf << "Deltas: " << DeltaS << " exp(-DeltaS) : "<< w << " Local update accepted at point: " << x << y << z << t << "\n";
}else{
//  op ++;
  return 0;
//  phi_new[s]=0;
//logf << "Deltas: " << DeltaS << " exp(-DeltaS) : "<< w <<" Local update rejected at point: " << x << y << z << t << "\n";
}
}

return 0;

}


void phi_update(int s, double eps){

double r[3];
ranlxd(r,3);
for(int a=0; a< Ng; a++ ){
phi_new[s][a] = phi[s][a] + eps*(r[a]*2 - 1.);
//std::cout << phi_new[s][a] << "\n";
}

}

double phi_DeltaS(int s, double lambda,double kappa){

  double delta =0;
  double  phi2=0;
  double  phi2new=0;
        for(int a = 0; a < Ng; a++){
        phi2 +=phi[s][a]*phi[s][a];
        phi2new +=phi_new[s][a]*phi_new[s][a];
        }


delta += phi2new - phi2 ;   // 1/2 m (phi'^2 -phi^2)
delta += lambda*(phi2new*phi2new - phi2*phi2 - 2*(phi2new -phi2 ) );

// Calculation of hopping part

int sp[Nd];
int sm[Nd];

su2 T[3];

_su2_t1(T[0]);
_su2_t2(T[1]);
_su2_t3(T[2]);

for(int mu= 0; mu < Nd; mu++){
sp[mu] = neig[s][mu];
sm[mu] = neig[s][Nd+mu];
}

for(int mu=0; mu< Nd; mu++){
su2 Um,Um_dn;

_su2_x_represent(Um,U[s][mu]);
_su2_x_represent(Um_dn,U[sm[mu]][mu]);

for(int a=0; a< Ng; a++ ){
for(int b=0; b< Ng; b++ ){

double W, W_dn;

su2 S1, S2, S3;

_suNg_times_suNg_dagger(S1,T[b],Um);
_suNg_times_suNg(S2,Um,S1);
_suNg_times_suNg(S3,T[a],S2);
_suNg_trace_re(W,S3);

_suNg_times_suNg_dagger(S1,T[b],Um_dn);
_suNg_times_suNg(S2,Um_dn,S1);
_suNg_times_suNg(S3,T[a],S2);
_suNg_trace_re(W_dn,S3);

delta -= 2*kappa*(phi_new[s][a] - phi[s][a])*W*phi[sp[mu]][b];
delta -= 2*kappa*phi[sm[mu]][a]*W_dn*(phi_new[s][b] - phi[s][b]);


}

}
}

return delta;


}


int phi_MC_test(int s, double DeltaS){

double t[1];
ranlxd(t,1);
double p,w;

if(DeltaS<0 || DeltaS ==0 ){
  for (int a = 0; a < Ng; a++) {
    phi[s][a] = phi_new[s][a];
  }
return 1;
}
if(DeltaS>0){
w=exp(-DeltaS);
p = t[0];
if(p<w){
  for (int a = 0; a < Ng; a++) {
    phi[s][a] = phi_new[s][a];
  }

return 1;
}else{
  return 0;
}
}

return 0;

}

void reunitarize(){

  for (int s = 0; s < V; s++) {

 for( int mu =0 ; mu < Nd; mu++){

  double detU;

  detU = U[s][mu].x[0]*U[s][mu].x[0] + U[s][mu].x[1]*U[s][mu].x[1] + U[s][mu].x[2]*U[s][mu].x[2] + U[s][mu].x[3]*U[s][mu].x[3];

//printf("detU at point %d :%f\n", s, detU );

U[s][mu].x[0]=U[s][mu].x[0]/sqrt(detU);
U[s][mu].x[1]=U[s][mu].x[1]/sqrt(detU);
U[s][mu].x[2]=U[s][mu].x[2]/sqrt(detU);
U[s][mu].x[3]=U[s][mu].x[3]/sqrt(detU);


}

}

}

void global_rotation(){

  double r[4];
  ranlxd(r,4);
  double alpha;

    su2 T[3];
    _su2_t1(T[0]);
    _su2_t2(T[1]);
    _su2_t3(T[2]);

  double R[3][3];
  int delta[3][3]={
    {1,0,0},
    {0,1,0},
    {0,0,1}
  };

  int eps[3][3][3] =
  {
  {{0,0,0},{0,0,1},{0,-1,0}},
  {{0,0,-1},{0,0,0},{1,0,0}},
  {{0,1,0},{-1,0,0},{0,0,0}}
};
//double  id[3][3];

  double n[3];
  su2 G;
  su2 G_adj;
 // su2 Id;

//  _suNg_unit(G);


  alpha = 2.* M_PI* r[0];

  n[0] = r[1]*10.-5.;
  n[1] = r[2]*10.-5.;
  n[2] = r[3]*10.-5.;

  double n_mod = sqrt(n[0]*n[0]+ n[1]*n[1] + n[2]*n[2]);

  n[0] = n[0]/n_mod;
  n[1] = n[1]/n_mod;
  n[2] = n[2]/n_mod;

//printf("%f %f %f\n",n.c[0],n.c[1],n.c[2]);

//ExpX(alpha, &n,&G);


(G).c[0].re = cos(alpha/2);
(G).c[0].im = n[2]*sin(alpha/2);
(G).c[1].re = n[1]*sin(alpha/2);
(G).c[1].im = n[0]*sin(alpha/2);
(G).c[2].re = -n[1]*sin(alpha/2);
(G).c[2].im = n[0]*sin(alpha/2);
(G).c[3].re = cos(alpha/2);
(G).c[3].im = -n[2]*sin(alpha/2);


_suNg_dagger(G_adj,G);

//_suNg_times_suNg(Id,G_adj,G);
//printf("%f + i %f  %f +i %f\n",Id.c[0].re, Id.c[0].im, Id.c[1].re, Id.c[1].im );
//printf("%f + i %f  %f +i %f\n",Id.c[2].re, Id.c[2].im, Id.c[3].re, Id.c[3].im );



for(int a = 0 ; a < 3 ; a++){
    for(int b = 0; b < 3 ; b++){
      R[a][b] = n[a]*n[b] + (delta[a][b]- n[a]*n[b])*cos(alpha);
      for(int c = 0; c < 3 ; c++){
      R[a][b] -= eps[a][b][c]*n[c]*sin(alpha);
        }

      }

    }


for(int s =0; s<V;s++)
  {
 double phi_d[3] ={0,0,0};

    for(int a = 0; a < 3 ; a++){
        for(int b = 0; b < 3 ; b++){
          phi_d[a] += R[b][a]*phi[s][b];
        }
      }

    phi[s][0] =phi_d[0];
    phi[s][1] =phi_d[1];
    phi[s][2] =phi_d[2];


for (int mu = 0; mu < 4; mu++) {
      su2 S2,S3,S4;
      _su2_x_represent(S4,U[s][mu]);
      _suNg_times_suNg(S2,S4,G_adj);
      _suNg_times_suNg(S3,G,S2);
      _su2_x_project(U[s][mu],S3);
      }

  }
}


void random_Z2() {

double r[2];
ranlxd(r,2);


if (r[0]>0.5) {

for(int s=0;s<V;s++){
      phi[s][0] = - phi[s][0];
      phi[s][1] = - phi[s][1];
      phi[s][2] = - phi[s][2];
  }
}

if (r[1]>0.5){

  int t, x, y, z;
  t = Nt-1;
  for (x=0; x<Nx; x++){
    for (y=0; y<Ny; y++){
      for (z=0; z<Nz; z++){

                int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;

                U[s][3].x[0]= -U[s][3].x[0];
                U[s][3].x[1]= -U[s][3].x[1];
                U[s][3].x[2]= -U[s][3].x[2];
                U[s][3].x[3]= -U[s][3].x[3];

              }
           }
        }
    }

}

void cycle(int nU, int multihit, int nphi, int nsweeps,std::ofstream &logf){

int j;

for(int i =0; i< nsweeps;i++){
j++;

update_gauge(nU,multihit);
update_phi(nphi);

eta = accepted/op*100.;

logf << "Acceptance rate at sweep " << i << " : " << eta << " %  \n";

if(eta<30)(eps = eps*0.98 );
if(eta>70)(eps = eps*1.02 );

}


  if(j == 10){
reunitarize();
    j=0;
  }

global_rotation();
random_Z2();

}

void update_gauge(int nU, int multihit){

double DeltaS;

for(int i =0; i< nU;i++){

  for (int s = 0; s < V; s++) {

      for(int mu =0; mu < Nd; mu++){
      // Calculation of the staple for the gauge sector
      su2_x st = staple(s,mu);

        for(int l =0; l < multihit;l++){
        // Generation of the new field U_new = X U, X close to identity
        gauge_update(s,mu,eps);
        // Calculation of Delta S for the gauge sector
        DeltaS = gauge_DeltaS(s,mu,st,beta,kappa);
        //std::cout << DeltaS << " "<< exp(-DeltaS) << " " << l << std::endl;
        // Metropolis test
        accepted += gauge_MC_test(s,mu,DeltaS);
        op ++;
        //plaqf  << avr_plaquette(neig)  << " " << i << "\n";
       }
      }
    }
  }
}

void update_phi(int nphi){
double DeltaS;


for(int i =0; i< nphi;i++){
    for (int s = 0; s < V; s++) {

phi_update(s, eps);

// Calculation of DeltaS

  DeltaS = phi_DeltaS(s,lambda, kappa);

// Metropolis test

accepted += phi_MC_test(s,DeltaS);
op++;

    }
  }

}

void U_copy(su2_x **U_out, su2_x **U_in){

for (int s = 0; s < V ; s++){
  for(int mu =0; mu < Nd; mu++){
    U_out[s][mu] = U_in[s][mu];
    }
  }

}

void phi_copy(double **phi_out, double **phi_in){

for (int s = 0; s < V ; s++){
  for(int a =0; a < Ng; a++){
    phi_out[s][a] = phi_in[s][a];
    }
  }

}
