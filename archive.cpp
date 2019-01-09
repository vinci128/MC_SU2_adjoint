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

  input in;

void read_input(char *argv[]){

  std::ifstream inputf;
  inputf.open(argv[1]);


   if (inputf.is_open())
    {
  inputf >> in.N;
  inputf >> in.b;
  inputf >> in.k;
  inputf >> in.l;
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


for(int s=0; s < V; s++ ){
	for(int a=0; a < Ng;a++){
		double buff = phi[s][a];

		infile.read((char*)&buff,sizeof(double));

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
