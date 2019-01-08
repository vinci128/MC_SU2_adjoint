
/*

Library for the implementation of complex numbers and SU(2) elements

Copyright: Vincenzo Afferrante 2018

*/

#ifndef __su2_h__
#define __su2_h__


// We define U = x_0 id + i x_i sigma_i

typedef struct{
double x[4];
} su2_x;

typedef struct{
  complex c[4];
} su2;



#define _su2_unit(u) \
_complex_1(u.c[0]);\
_complex_0(u.c[1]);\
_complex_0(u.c[2]);\
_complex_1(u.c[3])

#define _su2_t1(u) \
_complex_0(u.c[0]);\
_complex_1(u.c[1]);\
u.c[1].re = 0.5*u.c[1].re;\
_complex_1(u.c[2]);\
u.c[2].re = 0.5*u.c[2].re;\
_complex_0(u.c[3])

#define _su2_t2(u) \
_complex_0(u.c[0]);\
_complex_minus_i(u.c[1]);\
u.c[1].im = 0.5*u.c[1].im;\
_complex_i(u.c[2]);\
u.c[2].im = 0.5*u.c[2].im;\
_complex_0(u.c[3])

#define _su2_t3(u) \
_complex_1(u.c[0]);\
u.c[0].re = 0.5*u.c[0].re;\
_complex_0(u.c[1]);\
_complex_0(u.c[2]);\
_complex_minus_1(u.c[3]);\
u.c[3].re = 0.5*u.c[3].re

#define _su2_i_sigma1(u) \
_complex_0(u.c[0]);\
_complex_i(u.c[1]);\
_complex_i(u.c[2]);\
_complex_0(u.c[3])

#define _su2_i_sigma2(u) \
_complex_0(u.c[0]);\
_complex_minus_1(u.c[1]);\
_complex_1(u.c[2]);\
_complex_0(u.c[3])

#define _su2_i_sigma3(u) \
_complex_i(u.c[0]);\
_complex_0(u.c[1]);\
_complex_0(u.c[2]);\
_complex_minus_i(u.c[3])

#define _su2_i_sigma1_adj(u) \
_complex_0(u.c[0]);\
_complex_minus_i(u.c[1]);\
_complex_minus_i(u.c[2]);\
_complex_0(u.c[3])

#define _su2_i_sigma2_adj(u) \
_complex_0(u.c[0]);\
_complex_1(u.c[1]);\
_complex_minus_1(u.c[2]);\
_complex_0(u.c[3])

#define _su2_i_sigma3_adj(u) \
_complex_minus_i(u.c[0]);\
_complex_0(u.c[1]);\
_complex_0(u.c[2]);\
_complex_i(u.c[3])

#define _su2_null(u)\
u.x[0] = 0;\
u.x[1] = 0;\
u.x[2] = 0;\
u.x[3] = 0;

// u = v

#define _su2_eq(u,v) \
u.x[0] = v.x[0];\
u.x[1] = v.x[1];\
u.x[2] = v.x[2];\
u.x[3] = v.x[3];

#define _su2_adj(u,v) \
u.x[0] = v.x[0];\
u.x[1] = - v.x[1];\
u.x[2] = - v.x[2];\
u.x[3] = - v.x[3];

// v = y z

#define _su2_times_su2(v,y,z) \
v.x[0] = y.x[0]*z.x[0] - y.x[1]*z.x[1] - y.x[2]*z.x[2] - y.x[3]*z.x[3];\
v.x[1] = y.x[0]*z.x[1] + y.x[1]*z.x[0] - y.x[2]*z.x[3] + y.x[3]*z.x[2];\
v.x[2] = y.x[0]*z.x[2] + y.x[2]*z.x[0] - y.x[3]*z.x[1] + y.x[1]*z.x[3];\
v.x[3] = y.x[0]*z.x[3] + y.x[3]*z.x[0] - y.x[1]*z.x[2] + y.x[2]*z.x[1];

// u += y+z

#define _su2_add(u,y,z)\
u.x[0] += y.x[0]+ z.x[0];\
u.x[1] += y.x[1]+ z.x[1];\
u.x[2] += y.x[2]+ z.x[2];\
u.x[3] += y.x[3]+ z.x[3];


// u = y +z

#define _su2_add_assign(u,y,z)\
u.x[0] = y.x[0]+ z.x[0];\
u.x[1] = y.x[1]+ z.x[1];\
u.x[2] = y.x[2]+ z.x[2];\
u.x[3] = y.x[3]+ z.x[3];

#define _su2_x_trace(r,u)\
r = 2*u.x[0];


// u = v_0 id + i v_i sigma_i

#define _su2_x_represent(u,v) \
(u).c[0].re = v.x[0] ;\
(u).c[0].im = v.x[3] ;\
(u).c[1].re = v.x[2] ;\
(u).c[1].im = v.x[1] ;\
(u).c[2].re = - v.x[2] ;\
(u).c[2].im = v.x[1] ;\
(u).c[3].re = v.x[0] ;\
(u).c[3].im = -v.x[3];

// u = v_0 id + i v_i sigma_i

#define _su2_x_project(u,v) \
(u).x[0] = 0.5*(v.c[0].re+ v.c[3].re) ;\
(u).x[1] = 0.5*(v.c[1].im+ v.c[2].im) ;\
(u).x[2] = 0.5*(v.c[1].re- v.c[2].re) ;\
(u).x[3] = 0.5*(v.c[0].im- v.c[3].im) ;

// u = r*v (r real)


#define _su2_x_mul_assign(u,r,v)\
u.x[0] = r*v.x[0];\
u.x[1] = r*v.x[1];\
u.x[2] = r*v.x[2];\
u.x[3] = r*v.x[3];

// u += r*v (r real)

#define _su2_x_mul_add(u,r,v)\
u.x[0] += r*v.x[0];\
u.x[1] += r*v.x[1];\
u.x[2] += r*v.x[2];\
u.x[3] += r*v.x[3];

// u = v_a t_a

#define _su2_alg_represent(u,v) \
(u).c[0].re = 0.5*v[2] ;\
(u).c[0].im = 0 ;\
(u).c[1].re = 0.5*v[1] ;\
(u).c[1].im = -0.5*v[2] ;\
(u).c[2].re = 0.5*v[1] ;\
(u).c[2].im = 0.5*v[2] ;\
(u).c[3].re = -0.5*v[2] ;\
(u).c[3].im = 0;

/* u=v^dagger */
#define _suNg_dagger(u,v) \
   _complex_star((u).c[0],(v).c[0]); \
   _complex_star((u).c[1],(v).c[2]); \
   _complex_star((u).c[2],(v).c[1]); \
   _complex_star((u).c[3],(v).c[3])

/* u=v*w */
#define _suNg_times_suNg(u,v,w) \
      _complex_mul((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_assign((u).c[0],(v).c[1],(w).c[2]); \
      _complex_mul((u).c[1],(v).c[0],(w).c[1]);\
      _complex_mul_assign((u).c[1],(v).c[1],(w).c[3]); \
      _complex_mul((u).c[2],(v).c[2],(w).c[0]);\
      _complex_mul_assign((u).c[2],(v).c[3],(w).c[2]); \
      _complex_mul((u).c[3],(v).c[2],(w).c[1]);\
      _complex_mul_assign((u).c[3],(v).c[3],(w).c[3])

/* u=v*w^+ */
#define _suNg_times_suNg_dagger(u,v,w) \
      _complex_mul_star((u).c[0],(v).c[0],(w).c[0]);\
      _complex_mul_star_assign((u).c[0],(v).c[1],(w).c[1]); \
      _complex_mul_star((u).c[1],(v).c[0],(w).c[2]);\
      _complex_mul_star_assign((u).c[1],(v).c[1],(w).c[3]); \
      _complex_mul_star((u).c[2],(v).c[2],(w).c[0]);\
      _complex_mul_star_assign((u).c[2],(v).c[3],(w).c[1]); \
      _complex_mul_star((u).c[3],(v).c[2],(w).c[2]);\
      _complex_mul_star_assign((u).c[3],(v).c[3],(w).c[3])

/* u=v^+*w */
#define _suNg_dagger_times_suNg(u,v,w) \
      _complex_mul_star((u).c[0],(w).c[0],(v).c[0]);\
      _complex_mul_star_assign((u).c[0],(w).c[2],(v).c[2]); \
      _complex_mul_star((u).c[1],(w).c[1],(v).c[0]);\
      _complex_mul_star_assign((u).c[1],(w).c[3],(v).c[2]); \
      _complex_mul_star((u).c[2],(w).c[0],(v).c[1]);\
      _complex_mul_star_assign((u).c[2],(w).c[2],(v).c[3]); \
      _complex_mul_star((u).c[3],(w).c[1],(v).c[1]);\
      _complex_mul_star_assign((u).c[3],(w).c[3],(v).c[3])

/* k=Re Tr (u) */
#define _suNg_trace_re(k,u) \
   (k)=_complex_re((u).c[0])+ \
       _complex_re((u).c[3])     

/* k=Im Tr (u) */
#define _suNg_trace_im(k,u) \
   (k)=_complex_im((u).c[0])+ \
       _complex_im((u).c[3])    
      
      
#endif

