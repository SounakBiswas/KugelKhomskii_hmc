#include "mkl.h"
#include "math.h"
#include "mt19937ar.h"
#define DTAU 0.1
#define U 1.0
#define LX 4
#define NSITES LX
#define BETA 1.0
#define T_FL 0.2
#define WUP 100
#define MCS 100
#define IFPRECON 1
#define zcsrsymv mkl_cspblas_zcsrsymv
#define zcsrgemv mkl_cspblas_zcsrgemv
#define dcomplex MKL_Complex16
//typedef struct {double re; double im;} dcomplex;

int M;
double dt;
double dtby2;
double twodt;
int nsites,nf,M;
int twonsites;
double beta;
double dtau;
double corr_fn;
double corr_fn_err;
int tsteps;
double cgsteps;
int cg_ctr;
int lx;
int call;
int ione;
long long int accept;
double givens_param[5];
//fields
//dcomplex *phi;
dcomplex *phiA;
dcomplex *phiB;
dcomplex *phiC;
//dcomplex *R;
dcomplex *RA;
dcomplex *RB;
dcomplex *RC;
double *p_mom;
double *x_hs;
double *dVdx;
double * Proj;
int *lambda;
//preconditioning
double *L_pcg;
double *D_pcg;
dcomplex *aux_pcg;
double alpha_pcg;
//* ******** */
double root2;
double onebyroot2;
double root2U;
//ns temp matrix
dcomplex *aux_ns;

//nf temp matrix
dcomplex *aux_mtmx_nf;
dcomplex *aux1_nf;
dcomplex *aux2_nf;
dcomplex *aux3_nf;

double *pcopy;
double *xcopy;

//O**(-1) Phi vectors
//dcomplex *X_old;
dcomplex *XA_old;
dcomplex *XB_old;
dcomplex *XC_old;
//X= O**(-1) phi
//dcomplex *X;
dcomplex *XA;
dcomplex *XB;
dcomplex *XC;


//Checkerboard decomposition of inverse Green function matrics
dcomplex *acsr_kxa;
dcomplex *acsr_kxb;
dcomplex *acsr_kxbp;
dcomplex *acsr_kxap;
dcomplex *acsr_basis;
int *rowIndex_kxa;
int *rowIndex_kxb;
int *rowIndex_basis;
int *cols_basis;
int *cols_kxa;
int *cols_kxb;

int site2bond[NSITES][2];
int bond2site[NSITES][2];
