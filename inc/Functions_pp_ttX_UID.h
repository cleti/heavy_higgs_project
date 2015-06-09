
/*! \file
  \brief Interface for the evaluation of the unintegrated Catani/Seymour dipoles for the processes gg-> ttbar + g and qg -> ttbar + q, c.f. arXiv:hep-ph/0201036.
*/ 

#ifndef FUNCTIONS_UID_H
#define FUNCTIONS_UID_H

#include <vector>
#include <utility>
#include <functional>


#include "Global.h"
#include "Flags.h"
#include "Makros.h"
#include "Functions_Shared.h"
#include "Functions_pp_ttX_V.h"
#include "PhaseSpace.h"
#include "HistArray.h"


double Eval_UID_GG(
		   const PS_2_3& ps,
		   HiggsModel& hm,
		   const ulong& flags,
		   const double& dist_norm = 0.0,
		   DistVec* dist = nullptr);

double Eval_UID_QG(
		   const PS_2_3& ps,
		   HiggsModel& hm,
		   const ulong& flags,
		   const double& dist_norm = 0.0,
		   DistVec* dist = nullptr);




#endif
