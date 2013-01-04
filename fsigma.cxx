#include <stdio.h>
#include <math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include<vector>

#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "growth.h"

#define e 2.71828
#define pi 3.141593
#define rhoc 2.77536627e11

double solve_deltac(double par , double par1, double par2);

int calc_mstar(double Omegam, double Omegab, double hubble, double ns, double sigma8, double w, double redshift, vector<double> &k1, vector<double> &Tk1, int kbinnum, double *Mstar, double *nuu, int numbin, double *M){


  double delta_c = solve_deltac(redshift, Omegam, w);
  cout << "deltac="<<delta_c<<endl;
  /* call growth here*/

  double Da, Dadz;
  double func[kbinnum], A, k[kbinnum], Tk[kbinnum];
  double rhom0= Omegam*rhoc*pow(1.0+redshift, 3);
   
 growth growth_fact(Omegam, Omegab, 0, w);
 growth_fact.growthfactor(redshift, &Da, &Dadz);

 cout<< "growth factor="<<Da<<endl;

 for(int j=0; j< kbinnum; j++){
   k[j]= k1[j]; Tk[j]= Tk1[j];
   //  cout<< k[j]<<" "<<Tk[j]<<endl;
 }
   /*calculate normalization of the power sectrum*/


        for(int j=0; j< kbinnum; j++)
          {
            func[j]= pow(k[j], ns)*Tk[j]*Tk[j]*k[j]*k[j]* pow(3*(sin(k[j]*8)-k[j]*8*cos(k[j]*8)) /pow((k[j]*8),3),2.0);

          }

        {
gsl_interp_accel *acc
  = gsl_interp_accel_alloc ();

 gsl_spline *spline
   = gsl_spline_alloc (gsl_interp_cspline,kbinnum);

 gsl_spline_init (spline, k, func, kbinnum);

  A =pow(sigma8,2)/(1.0/(2.0*pi*pi)* gsl_spline_eval_integ (spline,k[0],k[kbinnum-1], acc));
  gsl_spline_free (spline);
 gsl_interp_accel_free (acc);
        }

        printf("A=%le\n", A);

double mhigh= 5e16;
double mlow= 5e11;
int binnum= 30;
double binsize= log(mhigh/mlow)/(binnum-1);
double mass[binnum], R, sigm[binnum], sigmsq; 

for(int i=0; i<binnum; i++){
  mass[i]= exp(log(mlow)+i*binsize); 
  R = pow(3.0*mass[i]/(4.0*pi*rhom0/pow(1+redshift,3.0)),0.3333);


                // calculate sigmsq and logD(sigmsq)as a function of M

 for(int j=0; j<kbinnum; j++)
   func[j]= 1.0/(2.0*pi*pi)*A*pow(k[j],ns)*Tk[j]*Tk[j]*k[j]*k[j]*pow(3*(sin(k[j]*R)-k[j]*R*cos(k[j]*R)) /pow((k[j]*R),3),2.0);


 {
gsl_interp_accel *acc
  = gsl_interp_accel_alloc ();

 gsl_spline *spline
   = gsl_spline_alloc (gsl_interp_cspline,kbinnum);

 gsl_spline_init (spline, k, func, kbinnum);

 sigmsq = pow(Da,2.0)*gsl_spline_eval_integ (spline,k[0],k[kbinnum-1], acc);
 sigm[i] = pow(sigmsq,0.5)/delta_c;                                                
 //cout<< i<<" "<<mass[i]<<" "<<sigm[i]<<endl;
 gsl_spline_free (spline);
 gsl_interp_accel_free (acc);
 }

 } 
//for(int i=0; i<numbin;i++) cout<<i<<" "<<M[i]<<endl;
 /*interpolate M vs sigm to get delc/sigm for mass*/
    {
         gsl_interp_accel *acc
           = gsl_interp_accel_alloc ();
         gsl_spline *spline
           = gsl_spline_alloc (gsl_interp_cspline, binnum);

         gsl_spline_init (spline, mass, sigm, binnum);

      for(int i=0;i<17;i++){
        if (M[i]>1){ nuu[i]= 1.0/gsl_spline_eval (spline, M[i], acc);}
          else {nuu[i]=0;}
	//cout<< i<<" "<<M[i]<<" "<<nuu[i]<<endl;
        }
          gsl_spline_free (spline);
          gsl_interp_accel_free (acc);
    }


/*interpolate sigm vs M to get M* */
/* double  sigm_i[binnum], mass_i[binnum];
 for (int i=0; i< binnum; i++){
   sigm_i[binnum-1-i]= sigm[i];
   mass_i[binnum-1-i]= mass[i];
   }

 {
         gsl_interp_accel *acc
           = gsl_interp_accel_alloc ();
         gsl_spline *spline
           = gsl_spline_alloc (gsl_interp_cspline, binnum);

	   gsl_spline_init (spline, sigm_i, mass_i, binnum);
	 //*Mstar= gsl_spline_eval (spline, 1.0, acc);
	 *Mstar=-1;
	           gsl_spline_free (spline); 
          gsl_interp_accel_free (acc);
     }
  
return 0;
*/
    *Mstar=-1;
}
