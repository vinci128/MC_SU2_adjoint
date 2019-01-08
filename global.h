#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "ranlxd.h"
#include "complex.h"
#include "su2.h"

#ifndef __global_h__

#define __global_h__

#define Ng 3
#define Nd 4
#define Nx 2
#define Ny 2
#define Nz 2
#define Nt 2
#define V Nx*Ny*Nz*Nt
#define Vd V*Nd
#define Vx Vd*4



extern su2_x U[V][Nd];
extern su2_x U_smear[V][Nd];
extern su2_x U_new[V][Nd];
extern su2_x X[V][Nd];
extern double phi[V][Ng];
extern double phi_smear[V][Ng];
extern double phi_new[V][Ng];


#endif
