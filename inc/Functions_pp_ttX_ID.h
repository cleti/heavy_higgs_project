

/*! \file
  \brief Interface for the evaluation of the integrated Catani/Seymour dipoles, c.f. arXiv:hep-ph/0201036, for the processes gg->ttbar and qg->ttbar. 
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


//! Evaluates the integrated dipoles that live on the normal 2->2 phase space, i.e. delta terms and distribution end-point terms.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_ID_GG(
	       const PS_2_2& ps,
	       HiggsModel& hm,
	       const ulong& flags);

//! Evaluates the integrated dipoles that live on the boosted 2->2 phase space, i.e. the + distribution terms.
/*!
  \param ps 2->2 phase space (with boosted initial state)
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_ID_X_GG(
		 const PS_2_2& ps_x,
		 HiggsModel& hm,
		 const ulong& flags);


double Eval_ID_X_QG(
		    const PS_2_2& ps_x,
		    HiggsModel& hm,
		    const ulong& flags);

#endif
