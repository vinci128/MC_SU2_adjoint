#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "ranlxd.h"
#include "complex.h"
#include "global.h"
#include "su2.h"
#include "geometry.h"

void fillneigh( int neigh[V][8] ){

int xup, yup, zup, tup;
int xmin, ymin, zmin, tmin;

for(int t=0; t < Nt; t++ ){

  tup = t+1;
  tmin = t-1;
  if(tup == Nx)(tup =0);
  if(tmin == -1)(tmin = Nt -1);


  for(int z=0; z < Nz; z++ ){

    zup = z+1;
    zmin = z-1;
    if(zup == Nz)(zup =0);
    if(zmin == -1)(zmin = Nz -1);


    for(int y=0; y < Ny; y++ ){

      yup = y+1;
      ymin = y-1;
      if(yup == Ny)(yup =0);
      if(ymin == -1)(ymin = Ny -1);

for(int x=0; x < Nx; x++ ){

xup = x+1;
xmin = x-1;
if(xup == Nx)(xup =0);
if(xmin == -1)(xmin = Nx -1);


  int s = x + Nx*y + Nx*Ny*z + Nx*Ny*Nz*t;

  int sp1 = xup + Ny*y + Ny*Nz*z + Ny*Nz*Nt*t;
  int sp2 = x + Ny*yup + Ny*Nz*z + Ny*Nz*Nt*t;
  int sp3 = x + Ny*y + Ny*Nz*zup + Ny*Nz*Nt*t;
  int sp4 = x + Ny*y + Ny*Nz*z + Ny*Nz*Nt*tup;

  int sm1 = xmin + Ny*y + Ny*Nz*z + Ny*Nz*Nt*t;
  int sm2 = x + Ny*ymin + Ny*Nz*z + Ny*Nz*Nt*t;
  int sm3 = x + Ny*y + Ny*Nz*zmin + Ny*Nz*Nt*t;
  int sm4 = x + Ny*y + Ny*Nz*z + Ny*Nz*Nt*tmin;

neigh[s][0] = sp1;
neigh[s][1] = sp2;
neigh[s][2] = sp3;
neigh[s][3] = sp4;
neigh[s][4] = sm1;
neigh[s][5] = sm2;
neigh[s][6] = sm3;
neigh[s][7] = sm4;

//printf("Position site: %d, neighbours%d %d %d %d %d %d %d %d\n", s, neigh[s][0],neigh[s][1],neigh[s][2],neigh[s][3],neigh[s][4],neigh[s][5],neigh[s][6],neigh[s][7]);

}
}
}
}


}
