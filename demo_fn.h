     struct quadratic_params
       {
         double rad_ratio, del_ratio, c;
       };
     
     double quadratic (double x, void *params);
     double quadratic_deriv (double x, void *params);
     void quadratic_fdf (double x, void *params, 
                         double *y, double *dy);

    struct vcirc_params
       {
         double vcirc, c;
       };

     double vcircfn (double x, void *params);
     double vcirc_deriv (double x, void *params);
     void vcirc_fdf (double x, void *params,
                         double *y, double *dy);

