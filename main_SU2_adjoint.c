#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "ranlxd.h"
#include "complex.h"
#include "global.h"
#include "geometry.h"
#include "update.h"
#include "archive.h"
#include "su2.h"
#include "APE_smearing.h"
#include "observables.h"


 su2_x U[V][Nd];
 su2_x U_smear[V][Nd];
 su2_x U_new[V][Nd];
 su2_x X[V][Nd];
 double phi[V][Ng];
 double phi_smear[V][Ng];
 double phi_new[V][Ng];

// Variables for the calculation of the acceptance rate

extern double accepted ;
extern double op;
extern double eta;

// Parameters of the theory (S = - beta/2 (1 - tr U mu nu(x)) )

extern double lambda;
extern double kappa ;
extern double beta;

// Parameter for the evolution

extern double eps ;


int main(int argc, char *argv[]){

for (int j =0; j < argc; j++){
std::cout << argv[j] << std::endl;
}
// Array of neighbours
int neig[V][8];
// Opening of the files for the observables
std::ofstream plaqf;
plaqf.open ("plaq.dat");
std::ofstream phif;
phif.open ("phi.dat");
std::ofstream logf;
logf.open ("log.dat");

std::ifstream inputf;
inputf.open(argv[1]);

std::string line;

 if (inputf.is_open())
  {
    while ( getline (inputf,line) )
    {
	    std::cout << line << '\n';
    }
    inputf.close();
  }

std::ofstream outf;
outf.open(argv[2]);

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

fillneigh(neig);

// Initialize the fields

cold_start();


logf << "Field initialized \n";

//std::cout << U[1][3].x[0] << "\n";

logf << "Update start \n ";


cycle(5,multihit,1,therm,logf);

// We use an adaptive algorithm for the update length


logf << "Thermalization completed" << std::endl;

// Calculation of the length and of the average plaquette (observables)
su2_x U_old[V][Nd];
double phi_old[V][Ng];

for(int i=0;i< n_meas;i++){
cycle(5,multihit,1,n_decorr,logf);

U_copy(U_old,U);
phi_copy(phi_old,phi);

plaqf << "no_smear: "<< avr_plaquette(neig)  << "\n";
phif  << "no_smear: " << phi_sq()  <<"\n";

for(int k =0; k < n_smear; k++){

APE_smearing(U_smear,U, 0.55);
APE_smearing_scalar(phi_smear,phi,U);

plaqf << "sm_level:" << k << " " << avr_plaquette_smear(neig)  << "\n";

phif  << "sm_level:" << k << " " << phi_sq_smear()  <<"\n";
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



plaqf.close();
phif.close();
logf.close();

logf << "Exiting program \n";

}

