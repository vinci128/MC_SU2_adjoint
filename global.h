#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "ranlxd.h"
#include "complex.h"
#include "su2.h"

#ifndef __global_h__

#define __global_h__

class input{
public:

  int N;
  double b;
  double k;
  double l;
  int th;
  int meas;
  int decorr;
  int mhit;
  int sm;
  int run;
  double alpha;
  int start;
  int init;

void read(char *argv[]);

};

extern input in;


extern  int Ng;
extern  int Nd;
extern  int Nx;
extern  int Ny;
extern  int Nz;
extern  int Nt;
extern  int V;
extern  int Vd;
extern  int Vx;

extern input in;

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

// Global array of neighbours

 extern int **neig;

// Global fields pointer


extern su2_x **U;
extern su2_x **U_smear;
extern su2_x **U_new;
extern su2_x **X;
extern double **phi;
extern double **phi_smear;
extern double **phi_new;






#endif
