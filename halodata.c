#include <stdio.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#include <math.h>

double start =1.0e12;  //bin start
double end= 1.0e17;     //bin end
int numbin= 20;         
long N= 300000;

int main(int argc, char **argv)
{   
    
  FILE *fin, *fout;
  char *file1, *file2;
  double hubble= atof(argv[1]);//0.71;
  double Omegam=atof(argv[2]);//0.264;
  double rhoc= 2.772e11; //Msun.h^2/Mpc^3
  double box_L= atof(argv[3]);//256.0;//512.0;   //Mpc/h
  int part_num= atoi(argv[4]);//256; //cube of this
  double minpart= atoi(argv[5]);//50; //smallest halo
  int particle_limit= atoi(argv[6]);//50; //min halo to consider
   file1= argv[7];
   file2= argv[8];
  int i, i1;
  double max, min, m[N];
  double minmass_corr, maxmass_corr, binsize_low, binsize_high, logbinsize; 
  int Nmtotal, Nm[numbin], count;
  double  mass[numbin], msum[numbin];

 logbinsize= (log(end)- log( start))/numbin;
 double sim_vol= pow(box_L/hubble,3); //Mpc^3
 double mass_res= Omegam*rhoc*pow(box_L/part_num,3); //Msun/h
// fin= fopen("/usr/projects/cosmo/sumanb/Gadget_w/halo_finder/output/outfnl500.mass","r");
// fout= fopen("/usr/projects/cosmo/sumanb/Gadget_w/halo_finder/output/outfnl500.mf256","w");
  fin= fopen(file1,"r");
  fout= fopen(file2,"w");
 i=0;
 max=1.0;
 min= 1.0e17;
 while (!feof(fin))
   {
     fscanf(fin,"%le",&m[i]);
     m[i]= m[i]*mass_res;
     if(m[i] >= max) max= m[i];
     if(m[i] <= min) min= m[i];
//if (m[i]==0.0) printf("%d %le\n",i, m[i]);     
i++;
   }
 count= i-1;
 fclose(fin);
 minmass_corr= mass_res*particle_limit;//100*5e11;
 maxmass_corr= max;

  printf("%d %le %le %lf %le\n", count, minmass_corr, maxmass_corr, logbinsize, mass_res);

 Nmtotal = 0;

 binsize_low = start;

 binsize_high = start*exp(logbinsize); 
 //printf("%le %le\n", binsize_low, binsize_high);
 fprintf(fout,"#bin no. mass [Msun/h] counts dn/d log m[1/Mpc^3], err(dn/d log m)\n");
 for( i1=0; i1<numbin; i1++){
   Nm[i1]=0;
   msum[i1]= 0.0;
   if(minmass_corr <= binsize_high && minmass_corr >= binsize_low){

     logbinsize= (log(binsize_high)- log(minmass_corr));
      //printf("logbinsize=%lf\n", logbinsize);
   }
 if(maxmass_corr <= binsize_high && minmass_corr >= binsize_low){

   logbinsize= (log(maxmass_corr)- log(binsize_low));
   //printf("logbinsize=%lf\n", logbinsize);
   }

 for (i=0; i<count; i++) {
  
   if(m[i] >= minmass_corr){
     if(m[i] < binsize_high && m[i] >= binsize_low){
       Nm[i1]= Nm[i1]+1;
       msum[i1]= msum[i1] + m[i];
       Nmtotal= Nmtotal+1;

     }
   }

}
 mass[i1]= msum[i1]/Nm[i1]; 

 if(msum[i1]==0 && Nm[i1]== 0) {
   mass[i1]= 0.0;
 } 
 // printf("%le %d %lf %lf %lf\n", mass[i1], Nm[i1], massmatch[i1], binsize_low, binsize_high);

if(Nm[i1]!=0) fprintf(fout, "%d %le %d %le %le\n ",i1, mass[i1], Nm[i1], Nm[i1]/logbinsize/sim_vol,sqrt(Nm[i1])/logbinsize/sim_vol);

 logbinsize= (log(end)-log(start))/numbin;

 binsize_low = binsize_high;

 binsize_high = binsize_low*exp(logbinsize);

 }

 printf("total halos= %d\n", Nmtotal);

 fclose(fout);
 return 0;      
}
 
