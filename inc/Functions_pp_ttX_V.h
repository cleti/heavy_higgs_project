
/*! \file
  \brief Interface for the evaluation of born level amplitudes and virtual corrections to the process pp -> ttbar + X.
*/ 

#ifndef FUNCTIONS_V_H
#define FUNCTIONS_V_H

#include <gsl/gsl_sf_dilog.h>

#include "Global.h"
#include "ScalarIntegrals.h"
#include "Flags.h"
#include "Makros.h"
#include "Functions_Shared.h"
#include "HiggsModel.h"


//! Evaluates the squared born amplitudes for the process \f$ gg \rightarrow \phi \rightarrow t\bar{t} \f$ and the interference with the QCD background at the given phase space point. NO spin/colour averaging factor included.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \param EFF Switch between effective Higgs-gluon coupling and full one-loop form factor.
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_B(
	      const PS_2_2& ps,
	      HiggsModel& hm,
	      const ulong& flags,
	      unsigned EFF=1);


//! Evaluates the squared born amplitudes for the process \f$ q\bar{q}  \rightarrow t\bar{t} \f$. NO spin/colour averaging factor included.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \sa PhaseSpace.h, HiggsModel.h
*/
double Eval_B_QQ(
		 const PS_2_2& ps,
		 HiggsModel& hm);


//////////////////////////////////////////////////////
double Eval_B_QCDxQCD(
		      const PS_2_2& ps,
		      const AmplitudePrefactors& ap,
		      const HiggsPrefactors& hp);

double Eval_B_2PHIxQCD(
		       const PS_2_2& ps,
		       const AmplitudePrefactors& ap,
		       const HiggsPrefactors& hp);

double Eval_B_PHIxPHI(
		      const PS_2_2& ps,
		      const AmplitudePrefactors& ap,
		      const HiggsPrefactors& hp);

double Eval_B_PHIxPHI_withINT12(
				const PS_2_2& ps,
				const HiggsModel& hm);
//////////////////////////////////////////////////////


//! Evaluates the virtual corrections to the process \f$gg \rightarrow \phi \rightarrow t\bar{t} \f$  and the itnerference with the QCD background  at the given phase space point. This includes the interference of one-loop QCD amplitudes with LO resonant amplitudes as well as that of LO QCD amplitudes with one-loop resonant amplitudes. Uses the effective Higgs-gluon vertex. NO spin/colour averaging factor included.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_V(
	      const PS_2_2& ps,
	      HiggsModel& hm,
	      const ulong& flags);



#endif
