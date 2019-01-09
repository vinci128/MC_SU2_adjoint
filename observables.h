#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "ranlxd.h"
#include "complex.h"
#include "global.h"
#include "geometry.h"
#include "update.h"
#include "su2.h"

#ifndef __observables_h__
#define __observables_h__

double plaquette(int s, int mu,int nu);
double plaq_smear(int s, int mu,int nu);
double avr_plaquette();
double avr_plaquette_smear();
double phi_sq();
double phi_sq_smear();

// Scalar operators, see Lee and Shigemitsu

void SC1_t(double O[]);
void SC2_t(double O[]);
void SC3_t(double O[]);

// Vector operators

void B_t(double O[][3]);

#endif
