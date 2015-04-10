
#ifndef FUNCTIONS_ID_H
#define FUNCTIONS_ID_H

#include <gsl/gsl_sf_dilog.h>

#include "Global.h"
#include "ScalarIntegrals.h"
#include "Flags.h"
#include "Makros.h"
#include "Functions_Shared.h"
#include "Functions_pp_ttX_V.h"



// integrated dipoles: delta terms and distribution end-point terms
double Eval_ID(
	       const PS_2_2& ps,
	       HiggsModel& hm,
	       const ulong& flags);
// integrated dipoles: distribution terms
double Eval_ID_X(
		 const PS_2_2& ps_x,
		 HiggsModel& hm,
		 const ulong& flags);


#endif
