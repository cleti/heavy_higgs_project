
#include "../../../../inc/Functions_Shared.h"

double Eval_B_EFF_PHIxPHI (
  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"
  
  double t1;
  double t11;
  double t15;
  double t3;
  double t5;
  double t8;
  t1 = s * s * s;
  t3 = FA0 * FA0;
  t5 = FH0 * FH0;
  t8 = At * At;
  t11 = Bt * Bt;
  t15 = DenS2(p1+p2, mH, GammaH);
  return(0.8e1 * t1 * (0.4e1 * t3 + t5) * (t11  + t8 * (0.1e1 - 0.4e1 / s ) )  * PREF_B_PHIxPHI * t15);
}
