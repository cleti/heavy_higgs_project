
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
#include "HistArray.h"


double Eval_UID_GG(
		   const PS_2_3& ps,
		   HiggsModel& hm,
		   const ulong& flags,
		   const double& dist_norm = 0.0,
		   std::vector<HistArray*>* dist = nullptr);

double Eval_UID_QG(
		   const PS_2_3& ps,
		   HiggsModel& hm,
		   const ulong& flags,
		   const double& dist_norm = 0.0,
		   std::vector<HistArray*>* dist = nullptr);




#endif