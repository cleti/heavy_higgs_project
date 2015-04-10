
#ifndef FUNCTIONS_R_H
#define FUNCTIONS_R_H

#include "Global.h"
#include "Flags.h"
#include "Makros.h"
#include "Functions_Shared.h"

double Eval_R_GG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);

double Eval_R_QQ(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);

double Eval_R_QG(
		 const PS_2_3& ps,
		 HiggsModel& hm,
		 const ulong& flags);


#endif
