/*bin the concentration in mass bin, compute the median c*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

using namespace std;

double  calc_c_vcirc(double vcirc_seek);
double calc_c_vcirc_cubic(double vcirc_seek);

int bin_conc(int tot_halos, int numbin, int particle_limit, double start, double end, double mass_res, vector<double> &m,vector<double> & FOFm, vector<double>& concentration, vector<double> &err_conc, vector<double> &r_ratio, vector<double> &err_r_ratio, vector<double> &v_circmax, vector<double> &err_vcircmax, vector<double> &v_vir, double *mass, double *FOFmass, int *Nm, double *rmedian, double *rerr, double *v_circmedian, double *v_circerr, double *c_vcirc, double *c_vcirc_err, double *cmedian, double *cmedianerr, double *cerr, double *cvarerrmean, double *cskew, double *ckurt, vector<double> &virial_rad)
{
  // FILE *fhist= fopen("c-M.hist","w");


  int i;
  double  max=1.0;
  double maxv=1.0;
  double minv=5000.0;
  //double min= 1.0e17;
  for(i=0; i< tot_halos; i++){
    if(FOFm[i] >= max) max= FOFm[i];
    if(v_circmax[i]*v_vir[i] >=maxv) maxv=v_circmax[i]*v_vir[i];    
    if(v_circmax[i]*v_vir[i] <minv && m[i]> mass_res*particle_limit) minv=v_circmax[i]*v_vir[i];
  }

  double minmass_corr= mass_res*particle_limit;//100*5e11;
  double maxmass_corr= max;
  cout<<"min FOF mass> "<<particle_limit <<"= "<<  minmass_corr<<" max FOF mass= "<<maxmass_corr<<"mass res.="<<mass_res<<"vmax="<<maxv<<" "<<"vmin="<< minv<<endl;
  int Nmtotal = 0;
  // start= minv; end=maxv;

  double logbinsize= (log(end)- log( start))/numbin;

  double binsize_low = start;

  double binsize_high = start*exp(logbinsize);
  //printf("%le %le\n", binsize_low, binsize_high);
  int i1;
  double msum[numbin], msumc[numbin], rsum[numbin], rvar[numbin],  FOFmsum[numbin],  csum[numbin], cvar[numbin], csumerr[numbin], logcsum[numbin], cvarerr[numbin];
  double v_circvar[numbin], v_circsum[numbin];
 

  for( i1=0; i1<numbin; i1++){
    Nm[i1]=0;
    msum[i1]= 0.0;
    FOFmsum[i1]=0.0;
    rsum[i1]= 0.0;
    rvar[i1]= 0.0;
    v_circsum[i1]=0.0;
    v_circvar[i1]=0.0;
    csum[i1]= 0.0;
    logcsum[i1]=0.0;
    cvar[i1]= 0.0;
    cvarerr[i1]= 0.0;
    csumerr[i1]= 0.0;
    msumc[i1]= 0.0;
    c_vcirc[i1]= 0.0;

    if(minmass_corr <= binsize_high && minmass_corr >= binsize_low){

      logbinsize= (log(binsize_high)- log(minmass_corr));
      //printf("logbinsize=%lf\n", logbinsize);
    }
    if(maxmass_corr <= binsize_high && minmass_corr >= binsize_low){

      logbinsize= (log(maxmass_corr)- log(binsize_low));
      //printf("logbinsize=%lf\n", logbinsize);
    }
    // cout<<"tot_halos before bin="<<tot_halos<<endl;
    for (i=0; i<tot_halos; i++) {
      if(/*FOF*/m[i] >= minmass_corr){
	if(/*FOF*/m[i] < binsize_high && /*FOF*/m[i] >= binsize_low /*&& v_circmax[i]*v_vir[i]< binsize_high &&  v_circmax[i]*v_vir[i] >= binsize_low*/ && concentration[i]>1.0 && concentration[i]<10.0 && virial_rad[i] > 0.0 && v_circmax[i]>=1.0 /*&& r_ratio[i]>0.0 && v_circmax[i]*v_vir[i]>minv*/){
 
	  // cout<<m[i]<<" "<<concentration[i]<<endl;
	  Nm[i1]= Nm[i1]+1;
	  msum[i1]= msum[i1] + m[i];
	  FOFmsum[i1]= FOFmsum[i1] + FOFm[i]; 
	  rsum[i1]= rsum[i1] + r_ratio[i]*m[i];
	  rvar[i1]= rvar[i1] + err_r_ratio[i]*m[i];
	  v_circsum[i1] += v_circmax[i]*m[i];
	  v_circvar[i1] +=  err_vcircmax[i]*m[i]; 
	  //c_vcirc[i1] += calc_c_vcirc_cubic(v_circmax[i])*m[i];
	  //c_vcirc_err[i1] += 0.5*(calc_c_vcirc_cubic(v_circmax[i]+ err_vcircmax[i])-calc_c_vcirc_cubic(v_circmax[i]- err_vcircmax[i]))*m[i];
	  csum[i1]= csum[i1] + concentration[i]*m[i];
	  logcsum[i1]= logcsum[i1] + log10(concentration[i])*m[i];
	  cvar[i1]= cvar[i1] + pow(concentration[i],2)*m[i];
	  csumerr[i1]= csumerr[i1] + err_conc[i]*m[i];
	  cvarerr[i1] += concentration[i]*m[i]*err_conc[i];//concentration[i];
	  Nmtotal= Nmtotal+1;
	  /*call for histrogram*/
	  //if (i1==15) fprintf(fhist,"%lf\n", concentration[i]);
          //cout<<m[i]<<" "<<concentration[i]<<" "<<calc_c_vcirc_cubic(v_circmax[i])<<endl;  
	}
      }
    }
    if(Nm[i1]>1){
      mass[i1]= (binsize_high+binsize_low)/2.0;//msum[i1]/Nm[i1];
    FOFmass[i1] = FOFmsum[i1]/Nm[i1];
    rmedian[i1]= rsum[i1]/msum[i1];
    rerr[i1]= sqrt(pow(rvar[i1]/msum[i1],2)+ pow(rmedian[i1],2)/Nm[i1]); 
    //pow(rvar[i1]/msum[i1]- rmedian[i1]*rmedian[i1],0.5);
  double fv=1.0;
    v_circmedian[i1]= fv*v_circsum[i1]/msum[i1];
    v_circerr[i1]= sqrt(pow(v_circvar[i1]/msum[i1],2)+pow(v_circmedian[i1],2)/Nm[i1]);//- v_circmedian[i1]*v_circmedian[i1],0.5);
    cmedian[i1]= csum[i1]/msum[i1];
    cmedianerr[i1]= pow(pow(csumerr[i1]/msum[i1],2)+pow(cmedian[i1]/sqrt(Nm[i1]), 2), 0.5);
    cerr[i1]=  pow(cvar[i1]/msum[i1] - pow(cmedian[i1],2),0.5)/cmedian[i1];
    cvarerrmean[i1]= cerr[i1]/sqrt(Nm[i1]);//sqrt(2*pow(cvarerr[i1]/msum[i1]/cerr[i1],2)+ pow(cerr[i1],2)/Nm[i1]);
  
    c_vcirc[i1] = calc_c_vcirc_cubic(v_circmedian[i1]);//c_vcirc[i1]/msum[i1]; 
    // cout<<Nm[i1]<<" "<<v_circmedian[i1]<<" "<<v_circerr[i1]<<endl;
    c_vcirc_err[i1] = 0.5*(calc_c_vcirc_cubic(v_circmedian[i1]+ v_circerr[i1])-calc_c_vcirc_cubic(v_circmedian[i1] - v_circerr[i1]));
    //cout<< calc_c_vcirc_cubic(1.1)<<endl;
    //v_circmedian[i1]<<" "<<v_circerr[i1]<<" "<<c_vcirc[i1]<<" "<<c_vcirc_err[i1]<<endl;
    cskew[i1]=0.0;
    ckurt[i1]=0.0;
    for (i=0; i<tot_halos; i++) {
      if(FOFm[i] >= minmass_corr){
	if(FOFm[i] < binsize_high && FOFm[i] >= binsize_low && concentration[i]>1.0){
	  cskew[i1] += pow(concentration[i],3)*m[i];
	  ckurt[i1] += pow(concentration[i]-cmedian[i1],4)*m[i];
	}
      }
    }
    //ckurt[i1]= ckurt[i1]/msum[i1]-4.0*cskew[i1]/msum[i1]*cmedian[i1]+6.0*pow(cerr[i1]*cmedian[i1],2)*pow(cmedian[i1],2)-3.0*pow(cmedian[i1],4.0);
    cskew[i1]= pow((cskew[i1]/msum[i1]-pow(cmedian[i1],3))/pow(cerr[i1]*cmedian[i1],3),1.0);///cmedian[i1];
    ckurt[i1]= pow(ckurt[i1]/msum[i1]/pow(cerr[i1]*cmedian[i1],4),1.0)-3;
    }
    if(msum[i1]==0 && Nm[i1]< 2) {
      mass[i1]= 0.0;
      rmedian[i1]=0.0;
      rerr[i1]= 0.0;
      cmedian[i1]= 0.0;
      cmedianerr[i1]=0.0;
      cerr[i1]= 0.0;
      cvarerrmean[i1]= 0.0;
      v_circmedian[i1]=0.0;
      v_circerr[i1]=0.0;
      c_vcirc[i1]=0.0;
      c_vcirc_err[i1]=0.0;
      cskew[i1]=0.0;
      ckurt[i1]=0.0;
       }			 
    // if(i1==8)  cout<<binsize_low<<" "<<binsize_high<<endl;
    logbinsize= (log(end)-log(start))/numbin;
    binsize_low = binsize_high;
    binsize_high = binsize_low*exp(logbinsize);
    }
   //fclose(fhist);
    cout<< "total halos >" << particle_limit<<" particles=" <<Nmtotal<<endl;
    return 0;
 }


    
  double  calc_c_vcirc_cubic(double vcirc_seek)
      {
	// Find the two bins that the density should be between
	// Interpolate the radius matching the requested density
	//cout<<vcirc_seek<<endl;
	 if(vcirc_seek<1.003) return 2.5;
	if (vcirc_seek>1.5) return 23.0;
       
	int possiblebin, bin, out;
	double c_vcirc;
	double ch= 24.0;
	double cl= 2.2;
	int nbin=99;
	double c_vcirc_i[nbin], vcirc_i[nbin];

	for (bin=0; bin< nbin; bin++){
	  c_vcirc_i[bin]= cl + (ch-cl)/(nbin-1)*bin;
	  vcirc_i[bin]= sqrt(0.216*c_vcirc_i[bin]/(log(1.0+c_vcirc_i[bin])-c_vcirc_i[bin]/(1.0+c_vcirc_i[bin])));
	  //cout<< c_vcirc_i[bin]<<" "<<vcirc_i[bin]<<endl;
	}
	{	
       gsl_interp_accel *acc 
          = gsl_interp_accel_alloc (); 
         gsl_spline *spline
           = gsl_spline_alloc (gsl_interp_cspline, nbin);
         gsl_spline_init (spline, vcirc_i, c_vcirc_i, nbin);
         c_vcirc= gsl_spline_eval (spline, vcirc_seek, acc);
 
          gsl_spline_free (spline);
         gsl_interp_accel_free (acc);

        return c_vcirc;
	}
      }

  double  calc_c_vcirc(double vcirc_seek)
      {
	// Find the two bins that the density should be between
	// Interpolate the radius matching the requested density
	if (vcirc_seek<1.0) return 0;
	int possiblebin, bin, out;
	double c_virc;
	double ch= 50.0;
	double cl= 1.0;
	int nbin=50;
	double c_vcirc_i[nbin], vcirc_i[nbin];

	for (bin=0; bin< nbin-1; bin++){
	  c_vcirc_i[bin]= cl + (ch-cl)/(nbin-1)*bin;
	  vcirc_i[bin]= sqrt(0.216*c_vcirc_i[bin]/(log(1.0+c_vcirc_i[bin])-c_vcirc_i[bin]/(1.0+c_vcirc_i[bin])));
	}

	for (bin = 0; bin < nbin-1; bin++) {
	  //cout<<delta_i[bin]<<" "<<r_i[bin]<<endl;
	  if (vcirc_i[bin] < vcirc_seek && vcirc_i[bin+1] > vcirc_seek) { 
	    possiblebin= bin; out=0;
	  }
	  if (vcirc_i[bin]==vcirc_seek) {
	    possiblebin= bin; out=1;
	  }
	} 
	if(possiblebin>nbin-1 ) return 0.0;
	/*linear interpolation*/
	//cout<< "possiblebin="<<possiblebin<<endl;
   
	if(out==0)  c_virc= c_vcirc_i[possiblebin]+ (c_vcirc_i[possiblebin+1]-c_vcirc_i[possiblebin])/(vcirc_i[possiblebin+1]- vcirc_i[possiblebin])*(vcirc_seek- vcirc_i[possiblebin]); 
	if(out==1) c_virc= c_vcirc_i[possiblebin];


	return c_virc;
      }
