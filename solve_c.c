     #include <stdio.h>
     #include <gsl/gsl_errno.h>
     #include <gsl/gsl_math.h>
     #include <gsl/gsl_roots.h>
     
     #include "demo_fn.h"
     #include "demo_fn.c"

     double solve_c (double param, double param1)
     {
       int status;
       int iter = 0, max_iter = 1000;
       const gsl_root_fsolver_type *T;
       gsl_root_fsolver *s;
       double r = 4.0, r_expected = 7.0;
       double x_lo = 0.01, x_hi = 50.0;
       gsl_function F;
       struct quadratic_params params = {param, param1, 0.0};
     
       F.function = &quadratic;
       F.params = &params;
     
       T = gsl_root_fsolver_brent;
       s = gsl_root_fsolver_alloc (T);
       gsl_root_fsolver_set (s, &F, x_lo, x_hi);
     
       do
         {
           iter++;
           status = gsl_root_fsolver_iterate (s);
           r = gsl_root_fsolver_root (s);
           x_lo = gsl_root_fsolver_x_lower (s);
           x_hi = gsl_root_fsolver_x_upper (s);
           status = gsl_root_test_interval (x_lo, x_hi,
                                            0.0, 0.01);
         }
       while (status == GSL_CONTINUE && iter < max_iter);
     
       gsl_root_fsolver_free (s);
     
       return r;
     }
