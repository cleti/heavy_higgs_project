
#include "../../../../inc/Functions_Shared.h"

double Eval_B_PHIxPHI (
  PS_2_2 const& ps)
{
  
#include "EVAL_V_PS_REFS"
  
  double t11;
  double t15;
  double t3;
  double t5;
  double t8;
  t3 = std::norm(F_ggH_p);
  t5 = std::norm(F_ggH_s);
  t8 = At * At;
  t11 = Bt * Bt;
  t15 = DenS2(p1+p2, mH, GammaH);
  return( (t3 + t5) * 0.5*s * (t11  + t8 * (0.1e1 - 0.4e1 / s ) ) * PREF_B_PHIxPHI * t15);
}
