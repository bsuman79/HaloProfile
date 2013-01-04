/* 
 * MINPACK-1 Least Squares Fitting Library
 *
 * Original public domain version by B. Garbow, K. Hillstrom, J. More'
 *   (Argonne National Laboratory, MINPACK project, March 1980)
 * 
 * Tranlation to C Language by S. Moshier (moshier.net)
 * 
 * Enhancements and packaging by C. Markwardt
 *   (comparable to IDL fitting routine MPFIT
 *    see http://cow.physics.wisc.edu/~craigm/idl/idl.html)
 */

/* Test routines for mpfit library
   $Id: testmpfit.c,v 1.3 2006/01/19 03:28:15 craigm Exp $
*/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf.h>

#include "mpfit.h"

#define  pi 3.141593
#define e 2.71828

/* This is the private data structure which contains the data points
   and their uncertainties */
struct vars_struct {
  double *x;
  double *y;
  double *ey;
};

/* Simple routine to print the fit results */
void printresult(double *x, double *xact, mp_result *result) 
{
  int i;

  if ((x == 0) || (result == 0)) return;
  printf("  CHI-SQUARE = %f    (%d DOF)\n", 
	 result->bestnorm, result->nfunc-result->nfree);
  printf("        NPAR = %d\n", result->npar);
  printf("       NFREE = %d\n", result->nfree);
  printf("     NPEGGED = %d\n", result->npegged);
  printf("     NITER = %d\n", result->niter);
  printf("      NFEV = %d\n", result->nfev);
  printf("\n");
  if (xact) {
    for (i=0; i<result->npar; i++) {
      printf("  P[%d] = %f +/- %f     (ACTUAL %f)\n", 
	     i, x[i], result->xerror[i], xact[i]);
    }
  } else {
    for (i=0; i<result->npar; i++) {
      printf("  P[%d] = %f +/- %f\n", 
	     i, x[i], result->xerror[i]);
    }
  }
    
}

/* 
 * linear fit function
 *
 * m - number of data points
 * n - number of parameters (2)
 * p - array of fit parameters 
 * dy - array of residuals to be returned
 * vars - private data (struct vars_struct *)
 *
 * RETURNS: error code (0 = success)
 */
//double get_A(double p[0]){ return &p[0]};
//double get_c(double p[1]){ return &p[1]};

int linfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey; 
  double conc;

  x = v->x;
  y = v->y;
  ey = v->ey;

  // n =2;
  // m=20;
  for (i=0; i<m; i++) {
  
    conc = p[0] - log(x[i]*p[1])- 2.0*log(1.0 + x[i]*p[1]);


    dy[i] = (y[i] - conc)/ey[i];

    }

  return 0;
}

/* Test harness routine, which contains test data, invokes mpfit() */
int testlinfit(double *x, double *y, double *ey, int m, int n, double *p, double *perror)
{
 
  p[0] = 1.0; p[1]= 2.0;          /* Initial conditions */
  double pactual[] = {1.0,1.0};   /* Actual values used to make data */
  //  double perror[n];                   /* Returned parameter errors */      
  int i;
  struct vars_struct v;
  int status;
  mp_par pars[n];                        /* Parameter constraints */
  mp_result result;

  bzero(&result,sizeof(result));       /* Zero results structure */
  result.xerror = perror;

  bzero(pars, sizeof(pars));           /* Initialize constraint structure */
  //pars[2].fixed = 1;                   /* Fix parameter 1 */

  v.x = x;
  v.y = y;
  v.ey = ey;
 
 /* Call fitting function for 10 data points and 2 parameters */
  status = mpfit(linfunc, m, n, p, pars, 0, (void *) &v, &result);

  //printf("*** testlinfit status = %d\n", status);
  //printresult(p, pactual, &result);
  
  return 0;
}

/* Main function which drives the whole thing */
int main(int argc, char *argv[])
{
  int i;
  int m= 20;
  int   n = 2;
  double p[n], perror[n];
  double x[m], y[m], ey[m];
   FILE *data;
   data = fopen("data.dat","r");
    
   for (i=0; i<m ; i++){
	fscanf(data,"%lf %lf %lf",&x[i], &y[i], &ey[i]);
      }
  
   testlinfit(x, y, ey, m, n, p, perror);
   printf("%lf %lf %lf %lf\n", p[0], perror[0], p[1], perror[1]);
 
  exit(0);
}
