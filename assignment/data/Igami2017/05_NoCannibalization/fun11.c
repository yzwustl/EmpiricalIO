#include "mex.h"  /* Always include this */
#include <math.h>  /* For mathematical functions */
#include "stdio.h" /* For mathematical functions */

/* Declare math functions */
double exp(double exponent);
double log(double antilog);
double pow(double base, double exponent);
        
/* The gateway routine */
void mexFunction(int nlhs, mxArray *plhs[],     /* Output variables */
        int nrhs, const mxArray *prhs[])        /* Input variables  */
{
    double *z11;
    double z1, z2;
    double beta, phi, kappa_inc, delta, year;
   
    /* Fool-proof for # Matlab inputs. */
    if(nrhs != 7)
        mexErrMsgTxt("fun11: 7 input arguments required.");
    
    /* Read inputs. */
    z1 =        mxGetScalar(prhs[0]);   /* EV_old for staying Old firm. */
    z2 =        mxGetScalar(prhs[1]);   /* EV_both for adopting Old firm. */
    beta =      mxGetScalar(prhs[2]);
    phi =       mxGetScalar(prhs[3]);
    kappa_inc = mxGetScalar(prhs[4]);
    delta =     mxGetScalar(prhs[5]);
    year =      mxGetScalar(prhs[6]);
        
    /* Initialize output. */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    z11 = mxGetPr(plhs[0]);
    
    /* Calculate the (3rd term of the) value of being a Corporate Venture. */
    *z11 = log(exp(0) 
                + exp(phi + beta * z1) 
                + exp(phi + beta * z2 - kappa_inc * pow(delta,year-1)));

    return;    
}
