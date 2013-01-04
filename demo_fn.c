 double
     quadratic (double x, void *params)
     {
       struct quadratic_params *p 
         = (struct quadratic_params *) params;
     
       double rad_ratio = p->rad_ratio;
       double del_ratio = p->del_ratio;
       double c = p->c;
     
       return  log(1.0 + x*rad_ratio) - x*rad_ratio/(1+x*rad_ratio) - log(1.0 + x)*del_ratio*pow(rad_ratio,3) + x/(1+x)*del_ratio*pow(rad_ratio,3) +c;
     }
     
     double
     quadratic_deriv (double x, void *params)
     {
       struct quadratic_params *p 
         = (struct quadratic_params *) params;
     
       double rad_ratio = p->rad_ratio;
       double del_ratio = p->del_ratio;
       double c = p->c;
     
       return -pow(rad_ratio,3)*del_ratio*x/pow(1+x,2) + pow(rad_ratio,2)*x/pow(1+ rad_ratio*x,2);
     }
     
     void
     quadratic_fdf (double x, void *params, 
                    double *y, double *dy)
     {
       struct quadratic_params *p 
         = (struct quadratic_params *) params;
     
       double rad_ratio = p->rad_ratio;
       double del_ratio = p->del_ratio;
       double c = p->c;
     
       *y = log(1.0 + x*rad_ratio) - x*rad_ratio/(1+x*rad_ratio) - log(1.0 + x)*del_ratio*pow(rad_ratio,3) + x/(1+x)*del_ratio*pow(rad_ratio,3) +c; 
       *dy = -pow(rad_ratio,3)*del_ratio*x/pow(1+x,2) + pow(rad_ratio,2)*x/pow(1+ rad_ratio*x,2);
     }
