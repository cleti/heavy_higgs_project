
/*! \file
  \brief Integrand functions used by VEGAS (wrapped in the Integrator class) to compute cross sections for the process pp -> H + X.
  \sa Integrator.h
*/ 



#ifndef INTEGRANDS_PPHX_H
#define INTEGRANDS_PPHX_H

#include "Global.h"
#include "Functions_Shared.h"
#include "Functions_pp_HX.h"


//! Integrand function for 2->1 scattering processes that contribute to pp -> H + X.
/*!
  \param x integration variables [ x[0], x[1]: PDF convolution ]
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_2_1_pdf(double* x, size_t dim, void* arg);

//! Integrand function for 2->2 scattering processes that contribute to pp -> H + X.
/*!
  \param x integration variables [ x[0], x[1]: PDF convolution, x[2]: scattering angle ]
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_2_2_pdf(double* x, size_t dim, void* arg);



#endif
