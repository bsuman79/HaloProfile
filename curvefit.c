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
#define rhoc 2.77536627e11

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


int intfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey; 
  double conc;
  x = v->x;
  y = v->y;
  ey = v->ey;
  for (i=0; i<m; i++) {
  
   //conc= p[0]*pow(p[1],3)/(log(1.0 + p[1]) - p[1]/(p[1]+1.0))/(x[i]*p[1])/pow(1+p[1]*x[i],2);
   // printf("%d %lf %lf\n", i, x[i],y[i]);    
    conc= p[0]*(log(1.0+p[1]*x[i]) - p[1]*x[i]/(1.0+p[1]*x[i]))/(log(1.0 + p[1]) - p[1]/(p[1]+1.0))/(4./3.*pi*pow(x[i],3)) ;
    //double err= 1/sqrt(ey[i])*y[i];
    if(y[i]!=0) dy[i] = (y[i] - conc)/(ey[i]);
   if(y[i]==0) dy[i] = 0.0;
   //if (1.0/y[i]==0) return 0; 
     //exit(EXIT_FAILURE);}
    }

  return 0;
}

int difffunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey;
  double conc;
  x = v->x;
  y = v->y;
  ey = v->ey;
  for (i=1; i<m; i++) {

    conc= p[0]*(log(1.0+p[1]*x[i]) - p[1]*x[i]/(1.0+p[1]*x[i])- log(1.0+p[1]*x[i-1]) + p[1]*x[i-1]/(1.0+p[1]*x[i-1]))/(log(1.0 + p[1]) - p[1]/(p[1]+1.0))/(4./3.*pi*(pow(x[i],3)- pow(x[i-1],3)));

    if(y[i]> 0) dy[i] = (y[i] - conc)/ey[i];
    }

  return 0;
}

int oneparamfunc(int m, int n, double *p, double *dy, double **dvec, void *vars)
{
  int i;
  struct vars_struct *v = (struct vars_struct *) vars;
  double *x, *y, *ey;
  double conc;
  x = v->x;
  y = v->y;
  ey = v->ey;
  for (i=1; i<m; i++) {

    conc= 200.0/3.0*p[0]*pow(p[1],3)/(x[i]*p[1])/pow(1+x[i]*p[1],2)/(log(1.0 + p[1]) - p[1]/(p[1]+1.0));
//(log(1.0+p[1]*x[i]) - p[1]*x[i]/(1.0+p[1]*x[i]))/(log(1.0 + p[1]) - p[1]/(p[1]+1.0))/(4./3.*pi*pow(x[i],3));

    dy[i] = (y[i] - conc)/ey[i];

    }

  return 0;
}



/* Test harness routine, which contains test data, invokes mpfit() */
int testlinfit(double *x, double *y, double *ey, int m, int n, double *p, double *perror)
{
 
  p[0] = 1.0; p[1]= 2.0;          /* Initial conditions */
  double pactual[] = {1.0,4.0};   /* Actual values used to make data */
  //  double perror[n];                   /* Returned parameter errors */      
  int i;
  struct vars_struct v;
  int status;
  mp_par pars[n];                        /* Parameter constraints */
  mp_result result;

  memset(&result,0,sizeof(result));  /* Zero results structure */
  result.xerror = perror;

  memset(pars,0,sizeof(pars));    /* Initialize constraint structure */
/* How to put limits on a parameter.  In this case, parameter 0 (log norm) is
     limited to be between 1.0 and 8.0. parameter 1 (conc) between 1.0 and 10.0 */
  //pars[0].limited[0] = 0;  pars[0].limited[1] = 1;
 // pars[0].limits[0] = 0.1; pars[0].limits[1] = 10.0;
   //pars[0].fixed = 1;
    pars[1].limited[0] = 0;  pars[1].limited[1] = 1;
    pars[1].limits[0] = 2.0; pars[1].limits[1] = 10.0;

  v.x = x;
  v.y = y;
  v.ey = ey;
 
 /* Call fitting function for 10 data points and 2 parameters */
  status = mpfit(difffunc, m, n, p, pars, 0, (void *) &v, &result);

  //printf("*** testlinfit status = %d\n", status);
  //printresult(p, pactual, &result);
  
  return 0;
}

/* Main function which drives the whole thing */
/*int main(int argc, char *argv[])
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
}*/
