
/*! \file
  \brief Leading order top/antitop decay matrix element.
*/ 

#ifndef FUNCTIONS_TDECAY
#define FUNCTIONS_TDECAY

#include "Functions_Shared.h"
#include "PhaseSpace.h"
#include "Lorentz.h"
#include "HiggsModel.h"


/* namespace ParametersTDecay */
/* { */
/*   extern double mW; */
/*   extern double mW2; */
/*   extern double GaW; */
/*   extern double GaW2;   */
/*   extern double GF; */
/*   extern double gw2; */
/*   extern double gw4; */
/*   extern void Init(const double& V, const double mScale); */
/* } */




double Eval_t_blnu(
		   const PS_2_3& ps,
		   FV& S,
		   const HiggsModel& hm
		   );


#endif
