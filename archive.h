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

#ifndef __archive_h__
#define __archive_h__

void read_input(char *argv[]);
void alloc_fields();

void write_adjoint_field(char filename[]);
void write_gauge_field(char filename[]);

void read_adjoint_field(char filename[]);
void read_gauge_field(char filename[]);


#endif
