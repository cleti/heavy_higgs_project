

#include <ctime>
#include "../inc/pVEGAS.h"
#include "../inc/Makros.h"


double exact = 1.3932039296856768591842462603255;

double
g (double *k, size_t dim, void *params)
{
  //std::this_thread::sleep_for( std::chrono::milliseconds(1) ); 
  double A = 1.0 / (M_PI * M_PI * M_PI);
  return A / (1.0 - cos (k[0]) * cos (k[1]) * cos (k[2]));
}


void
display_results (char *title, double result, double error)
{
  printf ("%s ==================\n", title);
  printf ("result = % .6f\n", result);
  printf ("sigma  = % .6f\n", error);
  printf ("exact  = % .6f\n", exact);
  printf ("error  = % .6f = %.2g sigma\n", result - exact,
          fabs (result - exact) / error);
}



int main()
{
  using namespace pVEGAS;


  int num_threads = 2;//omp_get_num_threads();
  int this_thread = omp_get_thread_num();

  PRINT(num_threads);
  PRINT(this_thread);
  
  double res, err;

  double xl[3] = { 0, 0, 0 };
  double xu[3] = { M_PI, M_PI, M_PI };

  const gsl_rng_type *T;
  gsl_rng **rngs = new gsl_rng*[num_threads];

  gsl_monte_function G = { &g, 3, 0 };

  size_t calls = 10000000;

  gsl_rng_env_setup ();
  gsl_rng_default_seed = time(0);
  T = gsl_rng_default;
  for (int i=0;i<num_threads;++i)
    {
      rngs[i] = gsl_rng_alloc (T);
    }


  {
    gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (num_threads,3);

    gsl_monte_vegas_integrate_parallel (&G, xl, xu, 3, calls/5/10, rngs, s,
                               &res, &err);
    display_results ("vegas warm-up", res, err);

    printf ("converging...\n");

    do
      {
	time_t INT_TIME = time(0);
        gsl_monte_vegas_integrate_parallel (&G, xl, xu, 3, calls/5, rngs, s,
                                   &res, &err);
        printf ("result = % .6f sigma = % .6f "
                "chisq/dof = %.1f\n", res, err, s->chisq);
	PRINT(difftime(time(0),INT_TIME));
      }
    while (fabs (s->chisq - 1.0) > 0.5);

    display_results ("vegas final", res, err);

    gsl_monte_vegas_free (s);
  }

    for (int i=0;i<num_threads;++i)
    {
      gsl_rng_free (rngs[i]);
    }
    delete []rngs;




  return 1;
}
