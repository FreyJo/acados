/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifdef CODEGEN_PREFIX
  #define NAMESPACE_CONCAT(NS, ID) _NAMESPACE_CONCAT(NS, ID)
  #define _NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else /* CODEGEN_PREFIX */
  #define CASADI_PREFIX(ID) vde_hess_chain_nm2_ ## ID
#endif /* CODEGEN_PREFIX */

#include <math.h>

#ifndef real_t
#define real_t double
#endif /* real_t */

#define to_double(x) (double) x
#define to_int(x) (int) x
#define CASADI_CAST(x,y) (x) y
/* Pre-c99 compatibility */
#if __STDC_VERSION__ < 199901L
real_t CASADI_PREFIX(fmin)(real_t x, real_t y) { return x<y ? x : y;}
#define fmin(x,y) CASADI_PREFIX(fmin)(x,y)
real_t CASADI_PREFIX(fmax)(real_t x, real_t y) { return x>y ? x : y;}
#define fmax(x,y) CASADI_PREFIX(fmax)(x,y)
#endif

#define PRINTF printf
#ifndef CASADI_SYMBOL_EXPORT
#if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
#if defined(STATIC_LINKED)
#define CASADI_SYMBOL_EXPORT
#else /* defined(STATIC_LINKED) */
#define CASADI_SYMBOL_EXPORT __declspec(dllexport)
#endif /* defined(STATIC_LINKED) */
#elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
#define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
#else /* defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__) */
#define CASADI_SYMBOL_EXPORT
#endif /* defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__) */
#endif /* CASADI_SYMBOL_EXPORT */
real_t CASADI_PREFIX(sq)(real_t x) { return x*x;}
#define sq(x) CASADI_PREFIX(sq)(x)

real_t CASADI_PREFIX(sign)(real_t x) { return x<0 ? -1 : x>0 ? 1 : x;}
#define sign(x) CASADI_PREFIX(sign)(x)

static const int CASADI_PREFIX(s0)[10] = {6, 1, 0, 6, 0, 1, 2, 3, 4, 5};
#define s0 CASADI_PREFIX(s0)
static const int CASADI_PREFIX(s1)[45] = {6, 6, 0, 6, 12, 18, 24, 30, 36, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
#define s1 CASADI_PREFIX(s1)
static const int CASADI_PREFIX(s2)[24] = {6, 3, 0, 6, 12, 18, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5};
#define s2 CASADI_PREFIX(s2)
static const int CASADI_PREFIX(s3)[7] = {3, 1, 0, 3, 0, 1, 2};
#define s3 CASADI_PREFIX(s3)
static const int CASADI_PREFIX(s4)[13] = {9, 1, 0, 9, 0, 1, 2, 3, 4, 5, 6, 7, 8};
#define s4 CASADI_PREFIX(s4)
static const int CASADI_PREFIX(s5)[49] = {45, 1, 0, 45, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44};
#define s5 CASADI_PREFIX(s5)
/* vde_hess_chain_nm2 */
CASADI_SYMBOL_EXPORT int vde_hess_chain_nm2(const real_t** arg, real_t** res, int* iw, real_t* w, int mem) {
  real_t a0=0.;
  if (res[0]!=0) res[0][0]=a0;
  if (res[0]!=0) res[0][1]=a0;
  if (res[0]!=0) res[0][2]=a0;
  real_t a1=arg[3] ? arg[3][0] : 0;
  if (res[0]!=0) res[0][3]=a1;
  a1=arg[3] ? arg[3][1] : 0;
  if (res[0]!=0) res[0][4]=a1;
  a1=arg[3] ? arg[3][2] : 0;
  if (res[0]!=0) res[0][5]=a1;
  a1=arg[3] ? arg[3][3] : 0;
  if (res[0]!=0) res[0][6]=a1;
  a1=arg[3] ? arg[3][4] : 0;
  if (res[0]!=0) res[0][7]=a1;
  a1=arg[3] ? arg[3][5] : 0;
  if (res[0]!=0) res[0][8]=a1;
  if (res[1]!=0) res[1][0]=a0;
  if (res[1]!=0) res[1][1]=a0;
  if (res[1]!=0) res[1][2]=a0;
  if (res[1]!=0) res[1][3]=a0;
  if (res[1]!=0) res[1][4]=a0;
  if (res[1]!=0) res[1][5]=a0;
  if (res[1]!=0) res[1][6]=a0;
  if (res[1]!=0) res[1][7]=a0;
  if (res[1]!=0) res[1][8]=a0;
  if (res[1]!=0) res[1][9]=a0;
  if (res[1]!=0) res[1][10]=a0;
  if (res[1]!=0) res[1][11]=a0;
  if (res[1]!=0) res[1][12]=a0;
  if (res[1]!=0) res[1][13]=a0;
  if (res[1]!=0) res[1][14]=a0;
  if (res[1]!=0) res[1][15]=a0;
  if (res[1]!=0) res[1][16]=a0;
  if (res[1]!=0) res[1][17]=a0;
  if (res[1]!=0) res[1][18]=a0;
  if (res[1]!=0) res[1][19]=a0;
  if (res[1]!=0) res[1][20]=a0;
  if (res[1]!=0) res[1][21]=a0;
  if (res[1]!=0) res[1][22]=a0;
  if (res[1]!=0) res[1][23]=a0;
  if (res[1]!=0) res[1][24]=a0;
  if (res[1]!=0) res[1][25]=a0;
  if (res[1]!=0) res[1][26]=a0;
  if (res[1]!=0) res[1][27]=a0;
  if (res[1]!=0) res[1][28]=a0;
  if (res[1]!=0) res[1][29]=a0;
  if (res[1]!=0) res[1][30]=a0;
  if (res[1]!=0) res[1][31]=a0;
  if (res[1]!=0) res[1][32]=a0;
  if (res[1]!=0) res[1][33]=a0;
  if (res[1]!=0) res[1][34]=a0;
  if (res[1]!=0) res[1][35]=a0;
  if (res[1]!=0) res[1][36]=a0;
  if (res[1]!=0) res[1][37]=a0;
  if (res[1]!=0) res[1][38]=a0;
  if (res[1]!=0) res[1][39]=a0;
  if (res[1]!=0) res[1][40]=a0;
  if (res[1]!=0) res[1][41]=a0;
  if (res[1]!=0) res[1][42]=a0;
  if (res[1]!=0) res[1][43]=a0;
  if (res[1]!=0) res[1][44]=a0;
  return 0;
}

CASADI_SYMBOL_EXPORT void vde_hess_chain_nm2_incref(void) {
}

CASADI_SYMBOL_EXPORT void vde_hess_chain_nm2_decref(void) {
}

CASADI_SYMBOL_EXPORT int vde_hess_chain_nm2_n_in(void) { return 5;}

CASADI_SYMBOL_EXPORT int vde_hess_chain_nm2_n_out(void) { return 2;}

CASADI_SYMBOL_EXPORT const char* vde_hess_chain_nm2_name_in(int i){
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    case 4: return "i4";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* vde_hess_chain_nm2_name_out(int i){
  switch (i) {
    case 0: return "o0";
    case 1: return "o1";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_hess_chain_nm2_sparsity_in(int i) {
  switch (i) {
    case 0: return s0;
    case 1: return s1;
    case 2: return s2;
    case 3: return s0;
    case 4: return s3;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const int* vde_hess_chain_nm2_sparsity_out(int i) {
  switch (i) {
    case 0: return s4;
    case 1: return s5;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int vde_hess_chain_nm2_work(int *sz_arg, int* sz_res, int *sz_iw, int *sz_w) {
  if (sz_arg) *sz_arg = 5;
  if (sz_res) *sz_res = 2;
  if (sz_iw) *sz_iw = 0;
  if (sz_w) *sz_w = 2;
  return 0;
}


#ifdef __cplusplus
} /* extern "C" */
#endif
