

/*! \file
  \brief Integrand functions used by VEGAS (wrapped in the Integrator class) to compute cross sections for the processes pp -> tt + X and pp -> tt + X -> l+ l- + b b-bar + jets.
  \sa Integrator.h
*/ 

#ifndef INTEGRANDS_PPTT_H
#define INTEGRANDS_PPTT_H

#include "Global.h"
#include "Makros.h"
#include "Functions_Shared.h"
#include "Integrator.h"

#include "Functions_pp_ttX_V.h"
#include "Functions_pp_ttX_ID.h"
#include "Functions_pp_ttX_R.h"
#include "Functions_pp_ttX_UID.h"

#ifdef WITH_T_SPIN
#include "Functions_tDecay.h"
#endif

//! Just for testing.
/*!
  \param x integration variables
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_poly2(double* x, size_t dim, void* arg);

//! Integrand without PDFs to reproduce some of W.B.'s plots.
/*!
  \param x integration variables
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_2_2(double* x, size_t dim, void* arg);

//! Integrand function for Born/virtual contribution to pp -> tt + X.
/*!
  \param x integration variables
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_2_2_pdf_BV(double* x, size_t dim, void* arg);

//! Integrand function for integrated dipole contribution to pp -> tt + X.
/*!
  \param x integration variables
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_2_2_pdf_ID(double* x, size_t dim, void* arg);

//! Integrand function for real corrections to gg -> tt + X.
/*!
  \param x integration variables
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_2_3_pdf(double* x, size_t dim, void* arg);

//! Integrand function for real corrections to qg/qqbar -> tt + X.
/*!
  \param x integration variables
  \param dim Integral dimension
  \param arg additional parameters
  \sa Integrator.h
*/
double Integrand_2_3_qg_qq_pdf(double* x, size_t dim, void* arg);

#endif
