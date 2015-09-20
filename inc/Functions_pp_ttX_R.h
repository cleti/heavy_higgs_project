
/*! \file
  \brief Interface for the evaluation of the real corrections to the process pp -> ttbar +X.
*/ 

#ifndef FUNCTIONS_R_H
#define FUNCTIONS_R_H

#include "Global.h"
#include "Flags.h"
#include "Makros.h"
#include "Functions_Shared.h"
#include "PhaseSpace.h"

//! Evaluates the squared tree-level amplitudes for the process \f$ gg \rightarrow \phi \rightarrow t\bar{t} + g \f$ and the interferences with the QCD background at the given phase space point. Uses the effective Higgs-gluon vertex. NO spin/colour averaging factor included. Leaves the Higgs propagator in HiggsModel in the state evaluated at s=(k1+k2)^2.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_R_GG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);

//! Evaluates the squared tree-level amplitudes for the process \f$ q\bar{q} \rightarrow \phi \rightarrow t\bar{t} + g \f$ and the interferences with the QCD background at the given phase space point. Uses the effective Higgs-gluon vertex. NO spin/colour averaging factor included.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_R_QQ(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);
//! Does not call SetHiggsPrefactors() on HiggsModel
double _Eval_R_QQ(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);

//! Evaluates the squared tree-level amplitudes for the process \f$ qg \rightarrow \phi \rightarrow t\bar{t} + q \f$ and the interferences with the QCD background at the given phase space point. Uses the effective Higgs-gluon vertex. NO spin/colour averaging factor included.
/*!
  \param ps 2->2 phase space
  \param hm Higgs model parameters
  \param flags specify which subamplitudes get evaluated
  \sa Flags.h, PhaseSpace.h, HiggsModel.h
*/
double Eval_R_QG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);
//! Does not call SetHiggsPrefactors() on HiggsModel
double _Eval_R_QG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);

#endif
