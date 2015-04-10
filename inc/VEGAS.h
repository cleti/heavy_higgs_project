/* monte/gsl_monte_vegas.h
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Michael Booth
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* header for the gsl "vegas" routines.  Mike Booth, May 1998 */

#ifndef __GSL_MONTE_VEGAS_H__
#define __GSL_MONTE_VEGAS_H__

/* configuration headers */
#include "VEGAS_config.h"

/* standard headers */
#include <math.h>
#include <stdio.h>

/* gsl headers */
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte.h>

#undef __BEGIN_DECLS
#undef __END_DECLS
#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
# define __BEGIN_DECLS /* empty */
# define __END_DECLS /* empty */
#endif

__BEGIN_DECLS

enum {GSL_VEGAS_MODE_IMPORTANCE      = 1, 
      GSL_VEGAS_MODE_IMPORTANCE_ONLY = 0, 
      GSL_VEGAS_MODE_STRATIFIED      = -1};

typedef struct {
  /* control variables - not changed after initialization */
  /* regarding general algorithm functionality */
  double alpha;
  int mode;
  int verbose;
  unsigned int iterations;
  int stage;

  /* control variables - not changed after initialization */
  /* regarding grid size */
  size_t dim;      // integral dimension
  size_t bins_max; // maximal size of grid
  double vol;      // total integation volume
  double * delx;   // integration intervals, upper - lower integration boundary ( size = dim )

  /* current grid and histogram parameters  */
  unsigned int bins;  // size of grid, bins along each axis -> may change between iterations
  unsigned int boxes; // number of boxes along each axis -> may change between iterations
  double * xi;     // histogram bin boundaries (grid) ( size = (bins+1)*dim ) -> within each bin the random points are sampled uniformly
  // the makro COORD refers to the array xi
  double * d;      // distribution of cumulated integrand values ( size = bins*dim )
  // the makro VALUE refers to the array d
  double * xin;    // temp., used for refining grid (axiswise, thus size = bins+1)
  double * weight; // temp., used for refining grid (axiswise, thus size = bins  )

  /* current event/point */
  double * x;      // uniformly distributed dim-dimensional random point for the next evaluation ...
  // bin and box coordinates = dim-dimensional vector of integers
  int    * bin;    // which lies in these bin ( size = dim )
  int    * box;    // and in this box         ( size = dim )
  double bin_vol;  // relative volume of this bin (remember that VEGAS transforms to [0,1] integration intervals)
  double wgt;      // statistical weight of the current random point = bin_vol*vol
  
  /* scratch variables preserved between calls to vegas1/2/3  */
  double jac;
  double wtd_int_sum; // weighted sum of results of each iteration
  double sum_wgts;    // sum of weights of each iteration (weight = 1/sigma^2 !)
  double chi_sum;
  double chisq;

  double result;
  double sigma;

  unsigned int it_start;
  unsigned int it_num;
  unsigned int samples;
  unsigned int calls_per_box;

  FILE * ostream;

} gsl_monte_vegas_state;

int gsl_monte_vegas_integrate(gsl_monte_function * f, 
                              double xl[], double xu[], 
                              size_t dim, size_t calls,
                              gsl_rng * r,
                              gsl_monte_vegas_state *state,
                              double* result, double* abserr);

gsl_monte_vegas_state* gsl_monte_vegas_alloc(size_t dim);

int gsl_monte_vegas_init(gsl_monte_vegas_state* state);

void gsl_monte_vegas_free (gsl_monte_vegas_state* state);

__END_DECLS

#endif /* __GSL_MONTE_VEGAS_H__ */

