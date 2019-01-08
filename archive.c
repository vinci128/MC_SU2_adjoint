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