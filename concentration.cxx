#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include<algorithm>
#include<iterator>

#include <math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

int npar=2;
#define rhoc 2.77536627e11
#define  pi 3.141593
#define mpc 3.0857e22 // m
#define msun 1.99e30 //kg
#define grav_c 6.673e-11 //m^3/kg/s^2

using namespace std;

int testlinfit(double*, double*, double*, int, int, double*, double*);

//double radius_ratio(double delta2, double delta1, int num_bin, double *delta_i, double *r_i);

double  calculateRadius(int num_bin, double delta_seek, double *delta_i, double *r_i, int *possiblebin);

int calc_concentration(double rhoz, string filename, int tot_halos, int num_col, int num_bin, int num_file, int *num_halos, double mass_res, int verbose, double delta1, double delta2, double delta3, vector<double> &r_ratio, vector<double> &err_r_ratio, vector<double> &v_circmax, vector<double> &err_vcircmax, vector<double> &v_vir, vector<double> &m, vector<double> &FOFm, vector<double> &concentration, vector<double> &err_conc, vector<double> &profile_norm, vector<double> &err_profile_norm, double delta_v, vector<double> &virial_rad ){

   int i, j, k, jj;
   // string  line;
   double temp[num_col], dist;
   double radius[num_bin], rho[num_bin], overden[num_bin], overdenint[num_bin], overdendiff[num_bin], vr[num_bin], rr, m1[num_bin], err_vcirc[num_bin], center_separation; 
   double delta_i[num_bin], r_i[num_bin], r_ii[num_bin], delta_diff[num_bin], r1, err_delta[num_bin];
   double r_iii[num_bin];
   int bin[num_bin];
   double count[num_bin];
   double p[npar], perror[npar];
   vector <double> v_circ(num_bin);
   FILE *fin;
   ostringstream sname, sname1;
   ifstream infile, infile1;
   tot_halos=0;
   int kk=0;
   int bad_halo=0;
   int nonvirial_halo=0;
   double r_v;
   if(verbose==0) dist=10.0;
   if(verbose==1) dist=0.07;
 for (i=0; i< num_file; i++){
   // cout<<"file open:"<<i<<endl;
   sname << filename <<".sodproperties."<< i;
   sname1 << filename <<".sodprofile."<< i;
  // infile.open(sname.str().c_str());
   fin= fopen(sname.str().c_str(),"r");
   infile1.open(sname1.str().c_str());
   // cout<< num_halos[i]<<endl; 
 for(j=0; j< num_halos[i]; j++){    
   //for(jj=0; jj< num_col; jj++)  infile>> temp[jj];// cout<<temp[15]<<endl;}
    fscanf(fin,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&temp[0], &temp[1], &temp[2], &temp[3], &temp[4], &temp[5], &temp[6], &temp[7], &temp[8], &temp[9], &temp[10], &temp[11], &temp[12], &temp[13], &temp[14], &temp[15], &temp[16], &temp[17], &temp[18], &temp[19], &temp[20]);

    for(jj=0; jj< num_bin; jj++){
     infile1>> bin[jj]; 
     infile1>> count[jj];
     infile1>> m1[jj];
     infile1>> radius[jj];
     infile1>> rho[jj];
     infile1>> overden[jj];
     infile1>> overdendiff[jj];
     infile1>> overdenint[jj];
     infile1>> vr[jj];
    }

    /* for(jj=0; jj< num_bin; jj++){
    overden[jj] *= rhoc/rhoz;
    overdendiff[jj] *= rhoc/rhoz;
    overdenint[jj] *= rhoc/rhoz;
    }*/
   double max_value=0;
   int max_pos=0;
   double countbin=0;
   int possiblebin;

   v_circ[0]= sqrt(4.0/3.0*M_PI*overdenint[0]*rhoc)*radius[0]*sqrt(grav_c/mpc*msun)*1e-3;
   for(jj=0; jj< num_bin; jj++){
   delta_i[jj]= overdenint[num_bin-1-jj];
   r_i[jj]= radius[num_bin-1-jj];   
   if(jj>0) v_circ[jj]= sqrt(4.0/3.0*M_PI*overdenint[jj]*rhoc)*radius[jj]*sqrt(grav_c/mpc*msun)*1e-3;      
   if (v_circ[jj]>max_value) {max_value= v_circ[jj]; max_pos= jj; countbin+= count[jj];}   
 }

   /*calculate r_Delta at z here from the density profile*/
     // r1=calculateRadius(num_bin, 200, delta_i, r_i, &possiblebin);
   FOFm[kk]= temp[1]*mass_res*(1.0-1.0/pow(temp[1],0.65));
   m[kk]= (double) /*temp[7];*/ temp[5]*mass_res;
   r1=pow(m[kk]/(4./3.*M_PI)/200.0/rhoc, 1.0/3.0);
   // cout<<r1<<" "<<temp[6]<<" "<<calculateRadius(num_bin, 200, delta_i, r_i, &possiblebin)<<endl;
   v_vir[kk]= sqrt(m[kk]/r1*grav_c/mpc*msun)*1e-3;
   center_separation= sqrt(pow(temp[8]-temp[11],2) + pow(temp[9]-temp[12],2) + pow(temp[10]-temp[13],2))/r1;


 /*determine v_circ max here for each halo*/

  if(max_value>v_vir[kk]){ v_circmax[kk]= max_value/v_vir[kk];}
    else {v_circmax[kk]=1.0;}
   
   err_vcircmax[kk]= v_circmax[kk]/sqrt(countbin);
   //cout<< m[kk]<<" "<<v_vir<<" "<<v_circmax[kk]<<endl;
  /*find when radial vel becomes less than -50 km/s*/ 

   double vrint[num_bin];
   double avgradvel;
   double rvirial;
   double totcount=0;

   vrint[0]=vr[num_bin-1];
   for (jj=1;jj<num_bin;jj++)  vrint[jj] = vr[num_bin-jj-1]+vrint[jj-1];
  
   for (jj=1; jj<num_bin-1;jj=jj++){
     //avgradvel = vrint[jj]/(jj+1);
      //cout<<jj<<" "<<avgradvel<<" "<<radius[num_bin-1-jj]<<endl;
      if(vr[num_bin-jj-1]>0.1*vr[num_bin-1]){
        rvirial= radius[num_bin-1-jj+1];  
        jj=num_bin-1;
      }
   }
   //cout<<rvirial<<endl;
   /*radius ratios*/
   int possiblebin1;
   double rad1, rad2;
   int totcount1, totcount2;
   double riii[num_bin], rii[num_bin];

 if(temp[5]/temp[1] > 0.1 && overdenint[2]> delta1 && center_separation< dist){
   rad1= calculateRadius(num_bin, delta1, delta_i, r_i, &possiblebin);
   rad2= calculateRadius(num_bin, delta2, delta_i, r_i, &possiblebin1);
 }  else { rad1=0.0; rad2=1.0;}

  totcount1=0;
  totcount2=0;
 for(jj=0; jj< num_bin-possiblebin; jj++) totcount1 += count[jj];
 for(jj=0; jj< num_bin-possiblebin1; jj++) totcount2 += count[jj];
 r_ratio[kk]= rad1/rad2; 
 err_r_ratio[kk]= sqrt(1.0/totcount1+ 1.0/totcount2)*r_ratio[kk];
 
 //cout<< r_ratio[kk]<<" "<<err_r_ratio[kk]<<" "<<possiblebin<< " "<<possiblebin1<<" "<<totcount1<<" "<<totcount2<<" "<<" "<<max_pos<<" "<<countbin<<endl;

   /*conc fit call here*/

 if(temp[5]/temp[1] > 0.1 && overdenint[3]> max(delta2, delta3) && center_separation < dist){
   int index;
   r_v= r1*calculateRadius(num_bin,delta_v,delta_i,r_i, &index)/calculateRadius(num_bin,200.0,delta_i,r_i, &index);
   // cout<<r1/r_v<<endl;

  if(delta3!=200.0) {
      double ratio= calculateRadius(num_bin, delta3, delta_i, r_i, &possiblebin)/calculateRadius(num_bin,200.0,delta_i,r_i, &possiblebin);
      r1=r1*ratio;
      m[kk]=m[kk]*pow(ratio,3)*delta3/delta2;
      v_vir[kk]=sqrt(m[kk]/r1);
      }
  /* for(jj=0; jj< num_bin; jj++){
    r_ii[jj]= radius[jj]/r1;
    totcount= count[jj] + totcount;
    delta_diff[jj]= overdenint[jj]*rhoc/m[kk];  
    err_delta[jj]= delta_diff[jj]/pow(totcount,0.5);
   }*/
   
  for(jj=0; jj< num_bin; jj++){
     rii[jj]=radius[jj]/r1;
     if(jj==0)  riii[jj]= 0.5*(rii[jj]+0.1);
     if(jj>0) riii[jj]=  0.5*(rii[jj]+rii[jj-1]);
   }
   
  index=0;
  double totcount=0;
  for(jj=0; jj< num_bin; jj++){    
    if(radius[jj]> 0.1*r_v && radius[jj]<= 1.0*r_v){
      r_ii[index]= rii[jj];
      r_iii[index]=riii[index];
      totcount= count[jj];//+totcount;
      delta_diff[index]= overdendiff[jj]*rhoc/m[kk];//*4.0/3.0*pi*pow(r_ii[index],3);
      err_delta[index]= delta_diff[index]/pow(totcount,0.5);
      // cout<<r_ii[index]<<" "<<overdendiff[index]<<" "<<overdenint[index]<<endl;//" "<<err_delta[index]<<" "<<totcount<<endl;
     index++; 
    }
  }
  // if(i==106 && j==192) for(jj=0; jj< num_bin; jj++) cout<<rii[jj]<<" "<<overdendiff[jj]<<endl;
  //cout<<"index="<<r_v<<endl;
  testlinfit(r_ii, delta_diff, err_delta, index, npar, p, perror);
 }  else { p[0]=0.0; p[1]=0.0;}
   profile_norm[kk]= p[0];
   err_profile_norm[kk]= perror[0];
   concentration[kk]= p[1];
   err_conc[kk]= perror[1];
 
   //cout<<i<<" "<<j<<" "<< m[kk]<<" "<<concentration[kk]<<" "<<err_conc[kk]<<endl;
    
   if(FOFm[kk]>=5.5e12 && FOFm[kk]<6.12e12 && concentration[kk]>1.0 && v_circmax[kk]>=1.0){
     // cout<<m[kk]<<" "<<v_circmax[kk]<<endl;
      for(jj=0; jj< num_bin; jj++){ 
	/*if (radius[jj]>=0.05*r_v && radius[jj]<=1.0*r_v && overdendiff[jj]>0){*/
         
	cout<<jj<<" "<<radius[jj]<<" "<<r_ii[jj]<<" "<<overdendiff[jj]<<" "<<count[jj]<<endl ;//<<" "<<v_circ[jj]<<endl;  //<<sqrt(overdenint[jj]/200.0)*r_ii[jj]/r1<<endl;     
 //	}
             }
	}
   // vmax^2/v200^2=(G*4/3pi*rhoc*delta*r^2)/(G*4/3pi*rhoc*200*r200^2)
   // vmax/v200=sqrt(delta/200)*r/r1

//cout<< center_separation<<" "<<r_v/rvirial<<endl;
   /*  if(m[kk]>1e15  && m[kk]<5e15 && concentration[kk]>1.0){
    cout<< "# c="<<concentration[kk]<<endl;  
    for(jj=0; jj< num_bin; jj++) cout<< r_ii[jj]<<" "<< overdendiff[jj]<<" "<<overdendiff[jj]/sqrt(count[jj])<<endl;
    }*/
   virial_rad[kk]= rvirial/r_v;
   if(concentration[kk]< 1.0){  
     //cout<<m[kk]<< " "<<concentration[kk]<<endl;
     bad_halo++;
   }
  
   if(virial_rad[kk] <1.0 && concentration[kk]>1.0) {
     nonvirial_halo++;
        }
   //if(j==0) { 
    //  for(jj=0; jj< num_bin; jj++) cout<< r_ii[jj]<<" "<<overdendiff[jj]<<" "<<overdendiff[jj]/sqrt(count[jj])<<endl;
     // } 
   // cout<<m[kk]<<" "<<v_circmax[kk]<<" "<<max_pos<<" "<<countbin<<" "<<v_circmax[kk]/sqrt(countbin)<<endl;
    kk++; 
   }
   //infile.clear();
   infile1.clear();
   //infile.close();
   infile1.close();
    fclose(fin);
   sname.str("");
   sname1.str("");
 }
 tot_halos=kk;
 cout<<"# concentrations calculated for "<< tot_halos<<" halos"<<endl;
 cout<<"# bad halos ="<< bad_halo<<" halos"<<endl;
 cout<<"# non_virial fraction ="<< nonvirial_halo<<" halos"<<endl;

 m.resize(tot_halos); r_ratio.resize(tot_halos);
 err_r_ratio.resize(tot_halos); FOFm.resize(tot_halos); 
 profile_norm.resize(tot_halos); err_profile_norm.resize(tot_halos);
 concentration.resize(tot_halos); err_conc.resize(tot_halos); 
 v_circmax.resize(tot_halos); err_vcircmax.resize(tot_halos);
 virial_rad.resize(tot_halos); v_vir.resize(tot_halos);
 return 0;
}


double calculateRadius(int num_bin, double delta_seek, double *delta_i, double *r_i, int *possiblebin)
{
  // Find the two bins that the density should be between
  // Interpolate the radius matching the requested density
  int bin, out;
  double r_seek;

  for (bin = 0; bin < num_bin-1; bin++) {
    //cout<<delta_i[bin]<<" "<<r_i[bin]<<endl;
    if (delta_i[bin] < delta_seek && delta_i[bin+1] > delta_seek) { 
      *possiblebin= bin; out=0;
     }
    if (delta_i[bin]==delta_seek) {
      *possiblebin= bin; out=1;
     }
  } 

  if(*possiblebin>19 ) return 0.0;
  /*linear interpolation*/
   //cout<< "possiblebin="<<possiblebin<<endl;
   
  if(out==0)  r_seek= r_i[*possiblebin]+ (r_i[*possiblebin+1]-r_i[*possiblebin])/(delta_i[*possiblebin+1]- delta_i[*possiblebin])*(delta_seek- delta_i[*possiblebin]); 
  if(out==1) r_seek= r_i[*possiblebin];

  return r_seek;
}



