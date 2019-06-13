#include <stdio.h>
#include <math.h>
#include "ranlxd.h"
#include "complex.h"
#include "global.h"
#include "geometry.h"
#include "update.h"
#include "su2.h"
#include "observables.h"
#include "archive.h"

#include <iostream>
#include <fstream>      // std::ofstream

using namespace std;


su2_x **U;
su2_x **U_smear;
su2_x **U_new;
su2_x **X;
double **phi;
double **phi_smear;
double **phi_new;

//  input in;


  void input::read(char *argv[]){

    std::ifstream inputf;
    inputf.open(argv[1]);

     if (inputf.is_open())
      {
    inputf >> N;
    inputf >> b;
    inputf >> k;
    inputf >> l;
    inputf >> th;
    inputf >> meas;
    inputf >> decorr;
    inputf >> mhit;
    inputf >> sm;
    inputf >> run;
    inputf  >> alpha;
    inputf >> start;
    }

        inputf.close();


  }



void alloc_fields(){

   U = new su2_x*[V];
   U_smear = new su2_x*[V];
   U_new = new su2_x*[V];
   X = new su2_x*[V];

   phi= new double*[V];
   phi_smear = new double*[V];
   phi_new = new double*[V];

  for(int i = 0; i < V; ++i){
      U[i] = new su2_x[Nd];
      U_smear[i] = new su2_x[Nd];
      U_new[i] = new su2_x[Nd];
      X[i] = new su2_x[Nd];

      phi[i] = new double[Nd];
      phi_smear[i] = new double[Nd];
      phi_new[i] = new double[Nd];
}
}

void dealloc_fields(){

  for(int i = 0; i < V; ++i){
      delete[] U[i];
      delete[] U_smear[i];
      delete[] U_new[i];
      delete[] X[i];

      delete[] phi[i];
      delete[] phi_smear[i];
      delete[] phi_new[i];
}

   delete[]  U;
   delete[] U_smear;
   delete[] U_new;
   delete[] X;

   delete[] phi;
   delete[] phi_smear;
   delete[] phi_new;


}


void write_adjoint_field(char filename[]){

  ofstream outfile;
  outfile.open(filename, ios::binary | ios::out);


for(int s=0; s < V; s++ ){
	for(int a=0; a < Ng;a++){
		double buff = phi[s][a];

		outfile.write((char*)&buff,sizeof(double));

	}
}
  outfile.close();

}


void write_gauge_field(char filename[]){

  ofstream outfile;
  outfile.open(filename, ios::binary | ios::out);

for(int s=0; s < V; s++ ){
	for(int mu=0; mu < Nd;mu++){
		su2 Um;

		_su2_x_represent(Um,U[s][mu]);

		outfile.write((char*)&Um,sizeof(su2));

	}
}

  outfile.close();

}


void read_adjoint_field(char filename[]){

  ifstream infile;
  infile.open(filename, ios::binary | ios::in);
	double buff;

for(int s=0; s < V; s++ ){
	for(int a=0; a < Ng;a++){


		infile.read((char*)&buff,sizeof(double));
//    std::cout << "buff" << buff << '\n';
phi[s][a] = buff;
	}
}
  infile.close();

}


void read_gauge_field(char filename[]){

  ifstream infile;
  infile.open(filename, ios::binary | ios::in);

for(int s=0; s < V; s++ ){
	for(int mu=0; mu < Nd;mu++){
		su2 Um;


		infile.read((char*)&Um,sizeof(su2));

		_su2_x_project(U[s][mu],Um);

	}
}

  infile.close();

}

void read_B(char *argv[], int n_smear){

double B[n_smear][Nt][3];

std::ifstream InFile;
InFile.open(argv[2], std::ios::in | std::ios::binary);
    for (int j = 0; j < n_smear; j++) {
for (int t = 0; t < Nt; t++) {
  for (int k = 0; k < 3; k++) {

  InFile.read( (char*)&B[j][t][k], sizeof(double));
  std::cout << B[j][t][k] << " " ;
    }
    std::cout  << '\n';
  }
  std::cout  << '\n';
}
std::cout  << '\n';
  InFile.close();
}
