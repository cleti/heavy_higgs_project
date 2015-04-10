
#ifndef FUNCTIONS_GGHX_H
#define FUNCTIONS_GGHX_H

#include "Global.h"
#include "Functions_Shared.h"
#include "ScalarIntegrals.h"




double Eval_B(FV const& p1,FV const& p2);
double Eval_B_eff(FV const& p1,FV const& p2);
double Eval_V(FV const& p1,FV const& p2);
double Eval_ID(FV const& p1,FV const& p2, double const& x);
double Eval_ID_X(FV const& p1,FV const& p2, double const& x);
double Eval_R_gg(FV const& p1,FV const& p2,FV const& p3);
double Eval_UID (FV const& p1,FV const& p2,FV const& p3);
double Eval_R_qq(FV const& p1,FV const& p2,FV const& p3);
// integrated in dim.reg. over divergent phase space and divergent part removed
double Eval_sigma_R_qg_fin(double const& z);
double Eval_sigma_R_qq_fin(double const& z);






#endif
