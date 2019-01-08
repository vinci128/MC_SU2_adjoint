#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "ranlxd.h"
#include "complex.h"
#include "global.h"
#include "geometry.h"
#include "su2.h"

#ifndef __update_h__
#define __update_h__

void cold_start();
void hot_start();

su2_x staple(int s, int mu,int neig[V][8]);
su2_x staple_in(su2_x in[][Nd],int s, int mu,int neig[V][8]);

void gauge_update(int s,int mu,double eps);
double gauge_DeltaS(int s, int mu,su2_x St, double beta,double kappa,int neig[V][8]);
int gauge_MC_test(int s, int mu,double DeltaS);

void phi_update(int s, double eps);
double phi_DeltaS(int s, double lambda,double kappa, int neig[V][8]);
int phi_MC_test(int s, double deltaS);

void reunitarize();

void global_rotation();
void random_Z2();


void cycle(int nU, int multihit, int nphi, int nsweeps,std::ofstream &logf);
void update_gauge(int nU, int multihit);
void update_phi(int nphi);

void U_copy(su2_x U_out[][Nd], su2_x U_in[][Nd]);
void phi_copy(double phi_out[][Ng], double phi_in[][Ng]);

#endif
