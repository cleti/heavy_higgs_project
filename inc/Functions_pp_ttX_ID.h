

/*! \file
  \brief Interface for the evaluation of the integrated Catani/Seymour dipoles for the processes gg->ttbar+g and qg->ttbar + q, c.f. arXiv:hep-ph/0201036. 
*/ 

#ifndef FUNCTIONS_ID_H
#define FUNCTIONS_ID_H

#include <gsl/gsl_sf_dilog.h>

#include "Global.h"
#include "ScalarIntegrals.h"
#include "Flags.h"
#include "Makros.h"
#include "Functions_Shared.h"
#include "Functions_pp_ttX_V.h"
#include "PhaseSpace.h"


//! Evaluates the integrated dipoles that live on the normal 2->2 phase space, i.e. delta terms and distribution end-point terms.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_ID_GG(
	       const PS_2_2& ps,
	       const double& x,
	       HiggsModel& hm,
	       const ulong& flags);

//! Evaluates the continuum part of the integrated dipoles that live on the boosted 2->2 phase space.
/*!
  \param ps 2->2 phase space (with boosted initial state)
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \param res_gg the result for the integrated dipoles associated with initial g->gg splitting
  \param res_qg the result for the integrated dipoles associated with initial q->qg splitting  
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
int Eval_ID_X(
	      const PS_2_2& ps_x,
	      const double& x,
	      HiggsModel& hm,
	      const ulong& flags,
	      double& res_gg,
	      double& res_qg);


#endif
