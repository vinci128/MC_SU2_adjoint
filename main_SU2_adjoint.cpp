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

// Opening of the files for the log
std::ofstream logf;
char log_name[128];
sprintf(log_name,"log_files/log_run%d_%dx%dx%dx%db%fk%fl%f",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda );
logf.open (log_name);

// Random number generator initialization

rlxd_init(1, in.init);

// Initialize the fields
if(in.start ==0){
  cold_start();
  logf << "Cold start\n";
} else{
  checkpoint_start();
  logf << "Start from configuration " << in.start << std::endl;
}

logf << "Field initialized \n";


logf << "Update start \n ";


cycle(5,multihit,1,therm,logf);

// We use an adaptive algorithm for the update length
logf << "Thermalization completed" << std::endl;


// Generation of configurations
int n_final = in.start + n_meas;
for(int i=in.start;i<n_final; i++){
cycle(5,multihit,1,n_decorr,logf);

sprintf(gauge_name,"conf/run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );
sprintf(adjoint_name,"adj_conf/adjoint_run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );

logf << gauge_name << std::endl;
logf << adjoint_name << std::endl;

write_gauge_field(gauge_name);
write_adjoint_field(adjoint_name);

}


dealloc_fields();

logf << "Arrays deallocated \n";

logf << "Exiting program \n";
logf.close();

std::ofstream checkfile;

char checkname[64];
sprintf(checkname,"test_run%d_k%f.out", in.run,kappa);
checkfile.open(checkname);

checkfile << "ok \n";

checkfile.close();


}
