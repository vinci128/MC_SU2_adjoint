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


void APE_smearing(su2_x out[][Nd], su2_x in[][Nd], double alpha);

void APE_smearing_scalar(double out[][Ng], double in[][Ng], su2_x gauge[][Nd]) ;

#endif