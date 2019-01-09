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

//input in;

int main(int argc, char *argv[]){

for (int j =0; j < argc; j++){
std::cout << argv[j] << std::endl;
}
// Opening of the files for the observables
std::ofstream plaqf;
plaqf.open ("plaq.dat");
std::ofstream phif;
phif.open ("phi.dat");
std::ofstream logf;
logf.open ("log.dat");

read_input(argv);

    beta = in.b;
    kappa = in.k;
    lambda = in.l;
     Nx = in.N;
     Ny = in.N;
     Nz = in.N;
     Nt = in.N;

V=Nx*Ny*Nz*Nt;

alloc_fields();

fillneigh();

std::ofstream outf;
outf.open(argv[2],ios::out|ios:binary);

// Parameters that determine the evolution

int therm = 1000;
int n_meas = 1000;
int n_decorr = 20;
int multihit = 10;
int n_smear = 5;
//int meas_freq = 10;
//int save_freq = 1000;
int n_run =1;
// Filename for configuration files


char gauge_name[128];
char adjoint_name[128];


// Random number generator initialization

rlxd_init(1, 12435);

// Initialize the fields

cold_start();


logf << "Field initialized \n";

//std::cout << U[1][3].x[0] << "\n";

logf << "Update start \n ";


cycle(5,multihit,1,therm,logf);

// We use an adaptive algorithm for the update length


logf << "Thermalization completed" << std::endl;


su2_x **U_old = new su2_x*[V];
double **phi_old = new double*[V];

for(int i = 0; i < V; ++i){
    U_old[i] = new su2_x[Nd];
    phi_old[i] = new double[Ng];
}
// Calculation of the length and of the average plaquette (observables)
for(int i=0;i< n_meas;i++){
cycle(5,multihit,1,n_decorr,logf);

U_copy(U_old,U);
phi_copy(phi_old,phi);

plaqf << "no_smear: "<< avr_plaquette()  << "\n";
phif  << "no_smear: " << phi_sq()  <<"\n";

for(int k =0; k < n_smear; k++){

  double B[Nt][3];

APE_smearing(U_smear,U, 0.55);
APE_smearing_scalar(phi_smear,phi,U);

plaqf << "sm_level:" << k << " " << avr_plaquette_smear()  << "\n";

phif  << "sm_level:" << k << " " << phi_sq_smear()  <<"\n";

B_t(B);
outf.write( (char*)&B, sizeof(B));

U_copy(U,U_smear);
phi_copy(phi,phi_smear);

}
U_copy(U,U_old);

phi_copy(phi,phi_old);

sprintf(gauge_name,"conf/run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );
sprintf(adjoint_name,"adj_conf/adjoint_run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );

write_gauge_field(gauge_name);
write_adjoint_field(adjoint_name);

}

logf << "Exiting program \n";


plaqf.close();
phif.close();
logf.close();
outf.close();

}
