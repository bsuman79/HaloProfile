#include <iostream>
#include <fstream>
#include <sstream>
#include <new>
#include <vector>
#include<iomanip>
#include<string>

#include <stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include <math.h>

#include "conc.h"


using namespace std;

int calc_concentration(double rhoz, string filename, int tot_halos, int num_col, int num_bin, int num_file, int *num_halos, double mass_res, int verbose, double delta1, double delta2, double delta3, vector<double> &r_ratio, vector<double> &err_r_ratio, vector<double> &v_circmax, vector<double> &err_vcircmax, vector<double> &v_vir, vector<double> &m, vector<double> &FOFm, vector<double> &concentration, vector<double> &err_conc, vector<double> &profile_norm, vector<double> &err_profile_norm, double delta_v, vector<double> &virial_rad);

double  solve_c (double param, double param1);

int bin_conc(int tot_halos, int numbin, int particle_limit, double start, double end, double mass_res, vector<double> &m,vector<double> & FOFm, vector<double>& concentration, vector<double> &err_conc, vector<double> &r_ratio , vector<double> &err_r_ratio, vector<double> &v_circmax, vector<double> &err_vcircmax, vector<double> &v_vir, double *mass, double *FOFmass, int *Nm, double *rmedian, double *rerr, double *v_circmedian, double *v_circerr, double *c_vcirc, double *c_vcirc_err, double *cmedian, double *cmedianerr, double *cerr, double *cvarerrmean, double *cskew, double *ckurt, vector<double> &virial_rad);

double calc_mstar(double, double, double, double, double, double, double, vector<double> &, vector<double> &, int, double*, double*, int, double*);

int main(int argc, char *argv[])
{ 
  int numbin, num_file, part_num, particle_limit;
  long tot_halos;
  string  filename, line;
  FILE   *fTk, *modelparams;
  char *file2, filenameTk[500];
  double hubble, Omegam, ns, sigma8, Omegab, redshift, w0, wa, wt, box_L;
  double delta1, delta2, delta3, rhoz;
  int model_no, verbose;

  cout<< "read in parameters"<<endl;
  cin >> numbin;
  cin >> num_file;
  cin >> model_no;
  cin >> redshift;
  cin >> box_L;
  cin >> part_num;
  cin >> particle_limit;
  cin >> delta1;
  cin >> delta2; 
  cin >> delta3;
  cin >> verbose; 
  cin >> filename;
  cin >> filenameTk;
  /*******************************************/
  // particle_limit=2000; //unfix the parameters
  //verbose=0;
  //redshift=0;
     /*****************************************/
  if(verbose==0) cout<< "#concentration calculated for relaxed+unrelaxed halos"<<endl;
  if(verbose==1) cout<< "#concentration calculated for relaxed halos only "<<endl;
  int i, i1, count, j, jj;
  double jnk1, jnk2, jnk3;
  ostringstream snameout;
  ofstream outfile;
  if(verbose==0) snameout<<"cM."<<box_L<<".f."<<model_no<<".z"<<redshift;
  if(verbose==1) snameout<<"cM."<<box_L<<".r."<<model_no<<".z"<<redshift;
  cout<<filename<<endl;
  outfile.open(snameout.str().c_str());  

  fTk = fopen(filenameTk,"r");
  vector<double> k(5000), Tk(5000);
  i=0;

  while (!feof(fTk) )
         {
          
           fscanf(fTk,"%le %lf %lf %lf %le",&k[i], &jnk1, &jnk2, &jnk3, &Tk[i]);
	    //printf("%le %le\n", k[i], Tk[i]);
  i++;
    }

    int  kbinnum = i - 1;
    printf("\n kbinnum = %d\n", kbinnum);
    k.resize(kbinnum); Tk.resize(kbinnum);
    fclose(fTk);

 modelparams=fopen("../halo_find/model_paramv1.txt","r");
 double par[num_models][num_par];

 for ( i=0; i< num_models; i++)
     {
       fscanf(modelparams,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&par[i][0], &par[i][1],&par[i][2], &par[i][3], &par[i][4], &par[i][5],&par[i][6], &par[i][7], &par[i][8], &par[i][9],&par[i][10], &par[i][11],&par[i][12]);
}

 fclose(modelparams);

        sigma8=par[model_no][4];
        ns=par[model_no][3];
        Omegam=par[model_no][6];
        Omegab=par[model_no][2]/par[model_no][8]/par[model_no][8];
        w0=par[model_no][5];
        hubble=par[model_no][8];
        wa=0.0;
        wt=w0+wa*redshift/(1+redshift);
        box_L= box_L*hubble; //convert box_L from Mpc to Mpc/h

      
  
  double logbinsize= (log(end)- log( start))/numbin;
  double Mstar, nuu[numbin];
  double mass[numbin], FOFmass[numbin], rmedian[numbin], rerr[numbin], v_circmedian[numbin], v_circerr[numbin], solvec[numbin], solvecerr[numbin], cmedian[numbin], cmedianerr[numbin], cerr[numbin], cvarerrmean[numbin];
  double c_vcirc[numbin], c_vcirc_err[numbin];
  double cskew[numbin], ckurt[numbin];
  int Nm[numbin];
  int num_halos[num_file];
  cout<<"redshift= "<<redshift<<" model= "<< model_no<<endl;
   double x= Omegam*pow(1+redshift,3)/(Omegam*pow(1+redshift,3)+ 1.0-Omegam)-1.0;
   if(delta3==-1) delta3= 18.0*pow(M_PI,2) +82.0*x-39.0*x*x;   
   if(delta3==-2) delta3= 200.0*Omegam*pow(1+redshift,3)/(Omegam*pow(1+redshift,3)+ 1.0-Omegam);
   cout<<"you chose virial radius, overdensity= "<<delta3<<endl;
   double  delta_v= 18.0*pow(M_PI,2) +82.0*x-39.0*x*x;
    double mass_res= Omegam*rhoc*pow(box_L/part_num,3); //Msun/h
   rhoz= rhoc*(Omegam*pow(1+redshift,3)+ 1.0-Omegam);
   ostringstream sname, sname1;
   ifstream infile;
   printf("\n parameters: %lf %lf %lf %lf %lf %lf %le\n", sigma8, ns, Omegam, wt, hubble, box_L, mass_res);
   tot_halos=0;
   for(i=0; i< num_file; i++){
   
   sname << filename <<".sodproperties."<< i;
   
   infile.open(sname.str().c_str());
   int k1=0;
   while(!infile.eof()) {
      getline(infile, line);
         k1++;
   }
   infile.clear();
   infile.close();
   sname.str("");
   num_halos[i] =k1-1;
   tot_halos= tot_halos+ num_halos[i];
   //cout<< tot_halos<<" "<<i<<endl;
}
  cout<< "# total halos= "<< tot_halos<<endl;

  vector<double> r_ratio(tot_halos), err_r_ratio(tot_halos), v_circmax(tot_halos), m(tot_halos), FOFm(tot_halos), concentration(tot_halos), err_conc(tot_halos), profile_norm(tot_halos), err_profile_norm(tot_halos), err_vcircmax(tot_halos), virial_rad(tot_halos);
  vector<double> v_vir(tot_halos); 
  cout<<"# calculating concentrations"<<endl;
  
  calc_concentration(rhoz, filename,tot_halos, num_col, num_bin, num_file, num_halos, mass_res, verbose, delta1, delta2, delta3, r_ratio, err_r_ratio, v_circmax, err_vcircmax, v_vir, m, FOFm, concentration, err_conc, profile_norm, err_profile_norm, delta_v, virial_rad);
 
  cout<<"# now binning the concentration measures in mass"<<endl;
 
  bin_conc(tot_halos, numbin, particle_limit, start, end, mass_res, m, FOFm, concentration, err_conc, r_ratio, err_r_ratio, v_circmax, err_vcircmax, v_vir, mass, FOFmass, Nm, rmedian, rerr, v_circmedian, v_circerr, c_vcirc, c_vcirc_err, cmedian, cmedianerr, cerr, cvarerrmean, cskew, ckurt, virial_rad); 

  //cout<< "# calculate M_star"<<endl;

  calc_mstar(Omegam, Omegab, hubble, ns, sigma8, wt, redshift, k, Tk, kbinnum, &Mstar, nuu, numbin, mass);

  cout<< "# mstar= "<<Mstar<<endl;

  // sname << "conc.z" <<redshift<<"."<< model_no;
  // ofstream outfile(sname.str().c_str(), ios::out);

   //outfile<<"#bin no.  FOFmass SOmass [Msun/h] SOmass/mstar counts rmedian (Mpc/h) err conc_median err conc_variance"<<endl;
   outfile<<" #bin no.  FOFmass SOmass [Msun/h] nu conc_median err conc_variance"<<endl;

   for( i1=0; i1<numbin; i1++){
     //cout<< i1<<" "<<nuu[i1]<<endl;
   //  if(Nm[i1]!=0) cout<<i1<<" "<< FOFmass[i1]<<" "<< mass[i1]<<" "<< mass[i1]/Mstar<<" "<< Nm[i1]<<" "<< rmedian[i1]<<" "<< rerr[i1]<<" "<< cmedian[i1]<<" "<< cmedianerr[i1]<<" "<< cerr[i1]<<" "<<nuu[i1]<<endl;
     double rplus= rmedian[i1]+ 0.05*rerr[i1];
     double rminus= rmedian[i1]- 0.05*rerr[i1];
     //double c_rad= solve_c(rmedian[i1], delta1/delta2);
     //double craderr= 0.5*(solve_c(rplus, delta1/delta2)-solve_c(rminus, delta1/delta2));
     if(Nm[i1]>2) outfile<< i1<<" "<< FOFmass[i1]<<" "<<mass[i1]<<" "<<nuu[i1]<<" "<<Nm[i1]<<" "<<cmedian[i1]/*1.35*/<<" "<<cmedianerr[i1]<<endl;//" "<<*/<<cerr[i1]/0.9<<endl;
//" "<<rmedian[i1]<<" "<<rerr[i1]<<" "<<c_rad<<" "<<craderr<<" "<<v_circmedian[i1]<<" "<<v_circerr[i1]<<" "<<c_vcirc[i1]<<" "<<c_vcirc_err[i1]<<endl;
    //<<" "<<cvarerrmean[i1]/*<<" "<<cskew[i1]<<" "<<ckurt[i1]*/<<" 
     //z1 1.4 z2 1.25 z0 1.08
}
    outfile.clear();
    outfile.close();
    snameout.str("");

  return 0;
}
