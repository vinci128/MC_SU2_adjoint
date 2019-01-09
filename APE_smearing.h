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
#include "observables.h"

#define M_PI 3.14159265358979323846

#ifndef __APE_smaering_h__
#define __APE_smaering_h__


void APE_smearing(su2_x **out, su2_x **in, double alpha);

void APE_smearing_scalar(double **out, double **in, su2_x **gauge) ;

#endif
