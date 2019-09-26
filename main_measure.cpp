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

char O1minus_name[128];

//  char plaq_name[128];
//  char phi_name[128];

//double B[Nt][3];

double B_p[Nt][3];

double B2[Nt][3];
double B_2p[Nt][3];
double Bphi[Nt][3];



sprintf(O1minus_name,"O1minus_output_files/output_run%d_Nt%d_Nx%d_Ny%d_Nz%d_B%f_K%f_L%f.bin",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda );
//sprintf(O0plus_name,"O0plus_output_files/output_Nt%d_Nx%d_Ny%d_Nz%d_B%f_K%f_L%f.bin",Nt,Nx,Ny,Nz, beta, kappa,lambda );

std::ofstream O1minusf;
O1minusf.open(O1minus_name,std::ios::out|std::ios::binary);

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


su2_x **U_old = new su2_x*[V];
double **phi_old = new double*[V];


for(int i = 0; i < V; ++i){
  U_old[i] = new su2_x[Nd];
  phi_old[i] = new double[Ng];
}

logf << "Update start \n ";


cycle(5,multihit,1,therm,logf);

// We use an adaptive algorithm for the update length
logf << "Thermalization completed" << std::endl;


// Generation of configurations for measurements
int n_final = in.start + n_meas;
for(int i=in.start;i<n_final; i++){
cycle(5,multihit,1,n_decorr,logf);

U_copy(U_old,U);
phi_copy(phi_old,phi);

U_copy(U_smear,U);
phi_copy(phi_smear,phi);

std::cout << "plaq:" << "no_smear: "<< avr_plaquette()  << "\n";
std::cout << "phi:"  << "no_smear: " << phi_sq()  <<"\n";

B_z_t(B_p);
B2_z_t(B2);
B_2z_t(B_2p);
Bphi_z_t(Bphi);

O1minusf.write((char*)&B_p, sizeof(B_p));
O1minusf.write((char*)&B2, sizeof(B2));
O1minusf.write((char*)&Bphi, sizeof(Bphi));
O1minusf.write((char*)&B_2p, sizeof(B_2p));

for(int k =0; k < n_smear; k++){



APE_smearing(U_smear,U, alpha);
APE_smearing_scalar(phi_smear,phi,U);


B_z_t(B_p);
B2_z_t(B2);
B_2z_t(B_2p);
Bphi_z_t(Bphi);

O1minusf.write( (char*)&B_p, sizeof(B_p));
O1minusf.write((char*)&B2, sizeof(B2));
O1minusf.write((char*)&Bphi, sizeof(Bphi));
O1minusf.write((char*)&B_2p, sizeof(B_2p));

U_copy(U,U_smear);
phi_copy(phi,phi_smear);

}

U_copy(U,U_old);
phi_copy(phi,phi_old);


}

// free dynamically allocated memory
for( int i = 0 ; i < V ; i++ )
{
    delete[] U_old[i]; // delete array within matrix
    delete[] phi_old[i];
}
// delete actual matrix
delete[] U_old;
delete[] phi_old;

// Last configurations are saved

sprintf(gauge_name,"conf/run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );
sprintf(adjoint_name,"adj_conf/adjoint_run%d_%dx%dx%dx%db%fk%fl%fn%d",n_run,Nt,Nx,Ny,Nz, beta, kappa,lambda,i );

logf << gauge_name << std::endl;
logf << adjoint_name << std::endl;

write_gauge_field(gauge_name);
write_adjoint_field(adjoint_name);

dealloc_fields();

logf << "Arrays deallocated \n";

logf << "Exiting program \n";
logf.close();
O1minusf.close();

std::ofstream checkfile;

char checkname[64];
sprintf(checkname,"test_run%d_k%f.out", in.run,kappa);
checkfile.open(checkname);

checkfile << "ok \n";

checkfile.close();


}
