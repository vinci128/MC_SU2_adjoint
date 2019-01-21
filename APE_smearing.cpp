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




void APE_smearing(su2_x **out, su2_x **in, double alpha)
{
double det;
for(int s =0; s < V;s++){
  for(int mu =0; mu<3; mu++ ) {

su2_x W,W_dag;

//printf("%f %f %f %f \n", in[s][mu].x[0], in[s][mu].x[1],in[s][mu].x[2],in[s][mu].x[3] );

_su2_x_mul_assign(out[s][mu],alpha,in[s][mu]);
W = staple_in(in,s,mu);
_su2_adj(W_dag,W);
_su2_x_mul_add(out[s][mu],(1.-alpha)/6.,W_dag );


det = out[s][mu].x[0]*out[s][mu].x[0] + out[s][mu].x[1]*out[s][mu].x[1] + out[s][mu].x[2]*out[s][mu].x[2] + out[s][mu].x[3]*out[s][mu].x[3];

_su2_x_mul_assign(out[s][mu],1./sqrt(det),out[s][mu]);

//printf("%f %f %f %f \n", out[s][mu].x[0], out[s][mu].x[1],out[s][mu].x[2],out[s][mu].x[3] );

    }
    _su2_eq(out[s][3],in[s][3]);

  }
}

void APE_smearing_scalar(double **out, double **in, su2_x **gauge) {

  su2 T[3];
	_su2_t1(T[0]);
	_su2_t2(T[1]);
	_su2_t3(T[2]);

int sp[Nd],sm[Nd];



su2 Um[4];
su2 U_adj[4];
su2 U_min[4];
su2 U_min_adj[4];

double W, Wmin;

for(int s =0; s < V;s++){

  for(int a =0; a <Ng; a++){
    out[s][a] = in[s][a]/7.;
  }

  for(int mu= 0; mu < 3; mu++){
sp[mu] = neig[s][mu];
sm[mu] = neig[s][Nd+mu];
}




for (int mu = 0; mu < 3; mu++) {

  _su2_x_represent(Um[mu],gauge[s][mu]); // U_mu(x)^(n-1)

  _suNg_dagger(U_adj[mu], Um[mu]); //U_mu(x)^dag

  _su2_x_represent(U_min[mu],gauge[sm[mu]][mu]); // U_mu(x)


  _suNg_dagger(U_min_adj[mu], U_min[mu]); //U_mu(x)^dag

  for (int a = 0; a < 3; a++) {

	 for (int b = 0; b < 3; b++) {

su2 S1, S2,S3;

//Construction of V_mu^(a b)(x) = Tr(T_a U_mu(x) T_b U^\dag_mu(x) )

        _suNg_times_suNg(S1,T[b],U_adj[mu]);
        _suNg_times_suNg(S2,Um[mu],S1);
        _suNg_times_suNg(S3,T[a],S2);
        _suNg_trace_re(W,S3);


        _suNg_times_suNg(S1,T[a],U_min_adj[mu]);
        _suNg_times_suNg(S2,U_min[mu],S1);
        _suNg_times_suNg(S3,T[b],S2);
        _suNg_trace_re(Wmin,S3);

W= 2*W;
Wmin = 2*Wmin;

out[s][a] += W*in[sp[mu]][b]/7. + Wmin*in[sm[mu]][b]/7.;

        }
      }


    }

  }

}
