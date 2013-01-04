     #include <stdio.h>
     #include<math.h>
     #include <gsl/gsl_errno.h>
     #include <gsl/gsl_matrix.h>
     #include <gsl/gsl_odeiv.h>
     

     int func (double a, const double y[], double f[],
           void *params)
     {
       double *p = (double *)params;
       double Om0= p[0];
       double w= p[1];
       double Ol0= 1.0-Om0;

       f[0] = y[1];
       f[1] = - (3./a -3./2.*(Om0/a+ Ol0*(1+w)/pow(a,3.0*w+1.0))/(Om0+Ol0/pow(a,3.0*w)))*y[1]+4./3.*y[1]*y[1]/(1.+y[0])+3./2.*Om0*pow(a,-5.0)/(Om0/pow(a,3.0)+Ol0/pow(a,3.0+3.0*w))*y[0]*(y[0]+1.0);
       return GSL_SUCCESS;
     }
     
     int jac (double a, const double y[], double *dfdy, 
          double dfda[], void *params)
     {
       double *p = (double *)params;
       double Om0= p[0];
       double w=  p[1];
       double Ol0= 1.0-Om0;

       gsl_matrix_view dfdy_mat 
         = gsl_matrix_view_array (dfdy, 2, 2);
       gsl_matrix * m = &dfdy_mat.matrix; 
       gsl_matrix_set (m, 0, 0, 0.0);
       gsl_matrix_set (m, 0, 1, 1.0);
       gsl_matrix_set (m, 1, 0,-4./3.*y[1]*y[1]/pow(1.+y[0],2.0)+3./2.*Om0*pow(a,-5.0)/(Om0/pow(a,3.0)+Ol0*pow(a,-3.0-3.0*w))*(2.0*y[0]+1.0) );
       gsl_matrix_set (m, 1, 1, -(3./a -3./2.*(Om0/a+ Ol0*(1+w)*pow(a,-3.0*w-1.0))/(Om0+Ol0/pow(a,3.0*w))) + 8./3.*y[1]/(1.+y[0]));
       dfda[0] = 0.0;
       dfda[1] = 0.0;
       return GSL_SUCCESS;
     }

 int func_lin (double a, const double y[], double f[],
           void *params)
     {
       double *p = (double *)params;
       double Om0= p[0];
       double w=  p[1];
       double Ol0= 1.0-Om0;

       f[0] = y[1];
       f[1] = - (3./a -3.0/2.0*(Om0/a+ Ol0*(1+w)/pow(a,3.0*w+1.0))/(Om0+Ol0/pow(a,3.0*w)))*y[1] + 3./2.*Om0*pow(a,-5.0)/(Om0/pow(a,3.0)+Ol0/pow(a,3.0+3.0*w))*y[0];
       return GSL_SUCCESS;
     }
     
     int jac_lin (double a, const double y[], double *dfdy, 
          double dfda[], void *params)
     {
       double *p = (double *)params;
       double Om0= p[0];
       double w= p[1];
       double Ol0= 1.0-Om0;

       gsl_matrix_view dfdy_mat 
         = gsl_matrix_view_array (dfdy, 2, 2);
       gsl_matrix * m = &dfdy_mat.matrix; 
       gsl_matrix_set (m, 0, 0, 0.0);
       gsl_matrix_set (m, 0, 1, 1.0);
       gsl_matrix_set (m, 1, 0,3./2.*Om0*pow(a,-5.0)/(Om0/pow(a,3.0)+Ol0*pow(a,-3.0-3.0*w))*(2.0*y[0]+1.0) );
       gsl_matrix_set (m, 1, 1, -(3./a -3./2.*(Om0/a+ Ol0*(1.0+w)*pow(a,-3.0*w-1.0))/(Om0+Ol0/pow(a,3.0*w))));
       dfda[0] = 0.0;
       dfda[1] = 0.0;
       return GSL_SUCCESS;
     }
     
   double solve_deltac(double par , double par1, double par2)
     {
       const gsl_odeiv_step_type * T 
         = gsl_odeiv_step_rk8pd;
     
       gsl_odeiv_step * s 
         = gsl_odeiv_step_alloc (T, 2);
       gsl_odeiv_control * c 
         = gsl_odeiv_control_y_new (1e-6, 0.0);
       gsl_odeiv_evolve * e 
         = gsl_odeiv_evolve_alloc (2);
       double ainitial= 1.0e-4;
     
       double p[2] = {par1, par2};
       double afinal= 1.0/(1+par);

       /*solve nonlinear equation to get the y_i s.t. y(z=z_obs)> 10^5*/
       gsl_odeiv_system sys = {func, jac, 2, &p};
       double yi;
       int i=0;
       double h = 1e-6;
       double yinitial= 1e-4;
       double ystep= 1e-7;
       double ymax= 5e5;
       double ybreak= 5e5;
       while (i<1e4)
      {
       double a = ainitial, a1 = afinal;
        yi= ystep*i + yinitial;
       double y[2] = { yi, 0 };
	   while(a< a1)
         {
           int status = gsl_odeiv_evolve_apply (e, c, s,
                                                &sys, 
                                                &a, a1,
                                                &h, y);
	    if (status != GSL_SUCCESS)
	      break;         

        }
         if (y[0]> ymax) {
           //printf ("%d %lf %.5e %.5e %.5e\n",i, a, yi, y[0], y[1]);	
           yinitial=  yi- ystep;
           ystep= ystep/10.;
           ymax= ymax*10;
           // printf("%le %le %le\n", yinitial, ystep, ymax);                             
           //break;
           i=0;
          }
         if (y[0] > ybreak) break;
          i++;
       } 

      /*now solve the linear equation to get delta_c at z=z_obs*/
      gsl_odeiv_system sys_lin = {func_lin, jac_lin, 2, &p};

       double a = ainitial, a1 = afinal;
             
       double y[2] = { yi, 0 };
	   while(a< a1)
         {
           int status = gsl_odeiv_evolve_apply (e, c, s,
                                                &sys_lin, 
                                                &a, a1,
                                                &h, y);
	    if (status != GSL_SUCCESS)
	      break;         
	   
       }
	   // printf ("%lf %.5e %.5e %.5e\n", a, yi, y[0], y[1]);          
       
      
       gsl_odeiv_evolve_free (e);
       gsl_odeiv_control_free (c);
       gsl_odeiv_step_free (s);
      
       return y[0];
     }

	/*int main(void){

	  int num_models= 38;
       double par[38][13];
       int i, ii;
       FILE* modelparams=fopen("model_paramv1.txt","r");

         for ( ii=0; ii< num_models; ii++)
     {
       fscanf(modelparams,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&par[ii][0], &par[ii][1],&par[ii][2], &par[ii][3], &par[ii][4], &par[ii][5],&par[ii][6], &par[ii][7], &par[ii][8], &par[ii][9],&par[ii][10], &par[ii][11],&par[ii][12]);
       printf("%lf %lf\n",par[ii][5], par[ii][6]);
}

  fclose(modelparams);

     double z[3]= {0, 1.0, 2.0};

     for(i=0; i< num_models; i++)   printf("deltac=%lf %lf %lf\n", solve_deltac(z[2], par[i][6], par[i][5]),solve_deltac(z[1], par[i][6], par[i][5]), solve_deltac(z[0], par[i][6], par[i][5]) );

     return;
	}*/
