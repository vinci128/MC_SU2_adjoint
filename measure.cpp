#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "ranlxd.h"
#include "complex.h"
#include "global.h"
#include "geometry.h"
#include "update.h"
#include "archive.h"
#include "su2.h"
#include "APE_smearing.h"
#include "observables.h"


int Nd = 4;
int Ng = 3;
int V;
int Vd = V*Nd;
int Vx = Vd*4;
int Nx;
int Ny;
int Nz;
int Nt;

input in;

int main(int argc, char *argv[]){


  in.read(argv);

      beta = in.b;
      kappa = in.k;
      lambda = in.l;
       Nx = in.N;
       Ny = in.N;
       Nz = in.N;
       Nt = in.N;

       V=Nx*Ny*Nz*Nt;

       alloc_fields();

       cold_start();


       std::cout << "Field initialized \n";
       fillneigh();
       // Parameters that determine the evolution

      // int therm = in.th;
       int n_meas = in.meas;
       //int n_decorr = in.decorr;
       //int multihit = in.mhit;
       int n_smear = in.sm;
       int n_run = in.run;
       double alpha = in.alpha;


// Filenames for configuration files

char gauge_name[128];
char adjoint_name[128];

// Open files for observables



  char O1minus_name[128];
  char O0plus_name[128];

  char plaq_name[128];
  char phi_name[128];

  //double B[Nt][3];

  double B_p[Nt][3];

  double T1[Nt][3];
  double T2[Nt][3];
  double T3[Nt][3];

  double SC1[Nt];
  double SC2[Nt];
  double SC3[Nt];

  sprintf(O1minus_name,"O1minus_output_files/output_Nt%d_Nx%d_Ny%d_Nz%d_B%f_K%f_L%f.bin",Nt,Nx,Ny,Nz, beta, kappa,lambda );
  sprintf(O0plus_name,"O0plus_output_files/output_Nt%d_Nx%d_Ny%d_Nz%d_B%f_K%f_L%f.bin",Nt,Nx,Ny,Nz, beta, kappa,lambda );

  sprintf(plaq_name,"obs/plaq_Nt%d_Nx%d_Ny%d_Nz%d_B%f_K%f_L%f.dat",Nt,Nx,Ny,Nz, beta, kappa,lambda );
  sprintf(phi_name,"obs/phi_Nt%d_Nx%d_Ny%d_Nz%d_B%f_K%f_L%f.dat",Nt,Nx,Ny,Nz, beta, kappa,lambda );


  std::ofstream O1minusf;
  O1minusf.open(O1minus_name,std::ios::out|std::ios::binary);

  std::ofstream O0plusf;
  O0plusf.open(O0plus_name,std::ios::out|std::ios::binary);

  std::ofstream plaqf;
  plaqf.open (plaq_name);
  std::ofstream phif;
  phif.open (phi_name);

  su2_x **U_old = new su2_x*[V];
  double **phi_old = new double*[V];


for(int i = 0; i < V; ++i){
    U_old[i] = new su2_x[Nd];
    phi_old[i] = new double[Ng];
}
// Calculation of the length, average plaquette and spectroscopical observables
for(int i=0;i< n_meas;i++){



  sprintf(gauge_name,"conf/run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );
  sprintf(adjoint_name,"adj_conf/adjoint_run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );

  std::cout << gauge_name << '\n';
  std::cout << adjoint_name << '\n';

  read_gauge_field(gauge_name);
  read_adjoint_field(adjoint_name);

  printf("configuration read n: " );
  printf("%d\n",i );


U_copy(U_old,U);
phi_copy(phi_old,phi);

U_copy(U_smear,U);
phi_copy(phi_smear,phi);

plaqf << "no_smear: "<< avr_plaquette()  << "\n";
phif  << "no_smear: " << phi_sq()  <<"\n";

std::cout << "plaq:" << "no_smear: "<< avr_plaquette()  << "\n";
std::cout << "phi:"  << "no_smear: " << phi_sq()  <<"\n";



//B_t(B);
B_z_t(B_p);
T1m_t(T1);
T2m_t(T2);
T3m_t(T3);

SC1_t(SC1);
SC2_t(SC2);
SC3_t(SC3);

O1minusf.write( (char*)&B_p, sizeof(B_p));
O1minusf.write( (char*)&T1, sizeof(T1));
O1minusf.write( (char*)&T2, sizeof(T2));
O1minusf.write( (char*)&T3, sizeof(T3));

O0plusf.write((char*)&SC1,sizeof(SC1));
O0plusf.write((char*)&SC2,sizeof(SC2));
O0plusf.write((char*)&SC3,sizeof(SC3));




for(int k =0; k < n_smear; k++){



APE_smearing(U_smear,U, alpha);
APE_smearing_scalar(phi_smear,phi,U);

plaqf << "sm_level:" << k << " " << avr_plaquette_smear()  << "\n";
phif  << "sm_level:" << k << " " << phi_sq_smear()  <<"\n";

B_z_t(B_p);
T1m_t(T1);
T2m_t(T2);
T3m_t(T3);

SC1_t(SC1);
SC2_t(SC2);
SC3_t(SC3);

O1minusf.write( (char*)&B_p, sizeof(B_p));
O1minusf.write( (char*)&T1, sizeof(T1));
O1minusf.write( (char*)&T2, sizeof(T2));
O1minusf.write( (char*)&T3, sizeof(T3));


/*
for(int t = 0; t < Nt; t++) {
  for (int l = 0; l < 3; l++) {
    printf("n= %d sm= %d  B[%d][%d] = %f \n",i,k+1, t,l,B[t][l] );
//    printf("n= %d sm= %d T1[%d][%d] = %f \n",i,k+1, t,l,T1[t][l] );
//    printf("n= %d sm= %d T2[%d][%d] = %f \n",i,k+1, t,l,T2[t][l] );
//    printf("n= %d sm= %d T3[%d][%d] = %f \n",i,k+1, t,l,T3[t][l] );

  }

}
*/


O0plusf.write((char*)&SC1,sizeof(SC1));
O0plusf.write((char*)&SC2,sizeof(SC2));
O0plusf.write((char*)&SC3,sizeof(SC3));

/*

for(int t = 0; t < Nt; t++) {
  printf("n= %d sm= %d SC1[%d] = %f \n",i,k, t,SC1[t] );
  printf("n= %d sm= %d SC2[%d] = %f \n",i,k, t,SC2[t] );
  printf("n= %d sm= %d SC3[%d] = %f \n",i,k, t,SC3[t] );
  printf("\n" );
}
  printf("\n" );
*/

U_copy(U,U_smear);
phi_copy(phi,phi_smear);

}
U_copy(U,U_old);

phi_copy(phi,phi_old);

}

plaqf.close();
phif.close();

O1minusf.close();
O0plusf.close();

}
