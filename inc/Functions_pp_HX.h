
/*! \file
  \brief Interface for the evaluation of matrix elements for the process pp -> H+X. 
*/ 


#ifndef FUNCTIONS_GGHX_H
#define FUNCTIONS_GGHX_H

#include "Global.h"
#include "Functions_Shared.h"
#include "ScalarIntegrals.h"
#include "HiggsModel.h"




//! Evaluate the squared Born matrix element for gg -> H. NO spin/colour averaging factor included.
/*!
  \param p1 incoming gluon 4-momentum
  \param p2 incoming gluon 4-momentum
  \param hm Model parameters
  \param EFF Switch between effective Higgs-gluon coupling and full one-loop form factor.
  \sa HiggsModel.h
*/
double Eval_B(FV const& p1,FV const& p2, HiggsModel const& hm, bool EFF);

//! Evaluate the virtual corrections for gg -> H: one-loop graphs with effective Higgs-gluon coupling. NO spin/colour averaging factor included.
/*!
  \param p1 incoming gluon 4-momentum
  \param p2 incoming gluon 4-momentum
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_V(FV const& p1,FV const& p2, HiggsModel const& hm);

//! Evaluate the integrated dipoles for gg -> H. These are the terms proportional to delta(1-x), c.f. arXiv:hep-ph/0201036. Uses the effective Higgs-gluon coupling. NO spin/colour averaging factor included.
/*!
  \param p1 incoming gluon 4-momentum
  \param p2 incoming gluon 4-momentum
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_ID(FV const& p1,FV const& p2, HiggsModel const& hm, double const& x);

//! Evaluate the integrated dipoles for gg -> H. This is the continuum part evaluated at the boosted phase space point, c.f. arXiv:hep-ph/0201036. Uses the effective Higgs-gluon coupling. NO spin/colour averaging factor included.
/*!
  \param p1 incoming gluon 4-momentum
  \param p2 incoming gluon 4-momentum
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_ID_X(FV const& p1,FV const& p2, HiggsModel const& hm, double const& x);

//! Evaluate the squared tree level matrix element for the process  gg -> H + g. Uses the effective Higgs-gluon coupling. NO spin/colour averaging factor included.
/*!
  \param p1 incoming gluon 4-momentum
  \param p2 incoming gluon 4-momentum
  \param p3 outgoing gluon 4-momentum
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_R_gg(FV const& p1,FV const& p2,FV const& p3, HiggsModel const& hm);

//! Evaluate the unintegrated dipoles for the process  gg -> H + g, c.f. arXiv:hep-ph/0201036. Uses the effective Higgs-gluon coupling. NO spin/colour averaging factor included.
/*!
  \param p1 incoming gluon 4-momentum
  \param p2 incoming gluon 4-momentum
  \param p3 outgoing gluon 4-momentum
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_UID (FV const& p1,FV const& p2,FV const& p3, HiggsModel const& hm);

//! Evaluate the squared tree level matrix element for the process  qqbar -> H + g. Uses the effective Higgs-gluon coupling. NO spin/colour averaging factor included.
/*!
  \param p1 incoming gluon 4-momentum
  \param p2 incoming gluon 4-momentum
  \param p3 outgoing gluon 4-momentum
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_R_qq(FV const& p1,FV const& p2,FV const& p3, HiggsModel const& hm);

//! Evaluate the squared tree level matrix element for the process  qg -> H + q (from Dawson, Nucl. Phys. B 359 (1991) 283-300). Uses the effective Higgs-gluon coupling. Spin/colour averaging as well as phase space factors are included.
/*!
  \param z  = mH2/s , mH: Higgs mass, s: partonic c.m.e.
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_sigma_R_qg_fin(double const& z, HiggsModel const& hm);

//! Evaluate the squared tree level matrix element for the process  qqbar -> H + g (from Dawson, Nucl. Phys. B 359 (1991) 283-300). Uses the effective Higgs-gluon coupling. Spin/colour averaging as well as phase space factors are included.
/*!
  \param z   = mH2/s , mH: Higgs mass, s: partonic c.m.e.
  \param hm Model parameters
  \sa HiggsModel.h
*/
double Eval_sigma_R_qq_fin(double const& z, HiggsModel const& hm);






#endif
