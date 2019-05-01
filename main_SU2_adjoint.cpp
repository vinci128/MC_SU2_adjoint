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

for (int j =0; j < argc; j++){
std::cout << argv[j] << std::endl;
}
// Opening of the files for the log
std::ofstream logf;
logf.open ("log.dat");


// We read the parameters from the input file

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

fillneigh();


// Parameters that determine the evolution

int therm = in.th;
int n_meas = in.meas;
int n_decorr = in.decorr;
int multihit = in.mhit;
//int n_smear = in.sm;
int n_run = in.run;
//double alpha = in.alpha;
// Filenames for configuration files

char gauge_name[128];
char adjoint_name[128];

// Random number generator initialization

rlxd_init(1, 12435);

// Initialize the fields

cold_start();


logf << "Field initialized \n";


logf << "Update start \n ";


cycle(5,multihit,1,therm,logf);

// We use an adaptive algorithm for the update length
logf << "Thermalization completed" << std::endl;


// Generation of configurations
for(int i=0;i< n_meas;i++){
cycle(5,multihit,1,n_decorr,logf);

sprintf(gauge_name,"conf/run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );
sprintf(adjoint_name,"adj_conf/adjoint_run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );

logf << gauge_name << std::endl;
logf << adjoint_name << std::endl;

write_gauge_field(gauge_name);
write_adjoint_field(adjoint_name);

}

logf << "Exiting program \n";
logf.close();


}
