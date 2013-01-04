#include <stdio.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <vector>

#define  pi 3.141593
#define mpc 3.0857e22 // m
#define msun 1.99e30 //kg
#define grav_c 6.673e-11 //m^3/kg/s^2

using namespace std;


int massfunction(int tot_halos, int numbin, double start, double end, double sim_vol,double min, double max, vector<double> &m, vector<double> &ms, double *mass, double* mass_s, int *Nm, double *mf, double *err_mf) {
 
  double rhoc= 2.77536627e11; //Msun.h^2/Mpc^3
  double minmass_corr, maxmass_corr, binsize_low, binsize_high;
  int Nmtotal, i, i1;
  double  msum[numbin], mssum[numbin];
  double logbinsize= (log(max)- log(min))/numbin;
  minmass_corr=min;
  maxmass_corr= max;
  printf("Min mass=%le Msun/h, Max mass=%le Msun/h\n", minmass_corr, maxmass_corr);

 Nmtotal = 0;

 binsize_low = start;

 binsize_high = start*exp(logbinsize); 
 //printf("%le %le\n", binsize_low, binsize_high);
 for( i1=0; i1<numbin; i1++){
   Nm[i1]=0;
   msum[i1]= 0.0;
   mssum[i1]= 0.0;
   if(minmass_corr <= binsize_high && minmass_corr >= binsize_low){

     logbinsize= (log(binsize_high)- log(minmass_corr));
      //printf("logbinsize=%lf\n", logbinsize);
   }
 if(maxmass_corr <= binsize_high && minmass_corr >= binsize_low){

   logbinsize= (log(maxmass_corr)- log(binsize_low));
   //printf("logbinsize=%lf\n", logbinsize);
   }

 for (i=0; i<tot_halos; i++) {
  
   if(m[i] >= minmass_corr){
     if(m[i] < binsize_high && m[i] >= binsize_low){
       Nm[i1]= Nm[i1]+1;
       msum[i1]= msum[i1] + m[i];
       mssum[i1] +=ms[i];
       Nmtotal= Nmtotal+1;

     }
   }
}
 mass[i1]= msum[i1]/Nm[i1]; 
 mass_s[i1]= mssum[i1]/Nm[i1];
 mf[i1]= Nm[i1]/logbinsize/sim_vol;
 err_mf[i1]= sqrt(Nm[i1])/logbinsize/sim_vol;
 if(msum[i1]==0 && Nm[i1]== 0) {
   mass[i1]= 0.0;
   mf[i1]=0.0;
   err_mf[i1]= 0.0;
   mass_s[i1]=0;
 } 
//if(Nm[i1]!=0) printf("%d %le %d\n ",i1, mass[i1], Nm[i1]);//, Nm[i1]/logbinsize/sim_vol,sqrt(Nm[i1])/logbinsize/sim_vol);

 logbinsize= (log(end)-log(start))/numbin;

 binsize_low = binsize_high;

 binsize_high = binsize_low*exp(logbinsize);

 }

  printf(" # halos = %d\n", Nmtotal);
     //cout<<m[i]<<endl;

 return 0;      
}
 
