
#include "../../../../inc/Functions_Shared.h"

double Eval_B_EFF_2PHIxQCD (
			    PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"

  double t1;
  double t11;
  double t17;
  double t2;
  double t20;
  double t23;
  double t26;
  double t3;
  double t39;
  c_double t43;
  double t44;
  double t7;
  t1 = sp(k2, p1);
  t2 = At * t1;
  t3 = sp(k1, p1);
  t7 = sp(p1, p2);
  t11 = t3 * t3;
  t17 = t3 * t7;
  t20 = t7 * t7;
  t23 = t11 * t3;
  t26 = Bt * t1;
  t39 = -0.2e1 * At * FH0 * t11 + At * FH0 * t23 - 0.2e1 * Bt * FA0 * t23 - 0.2e1 * FA0 * t11 * t26 - 0.2e1 * FA0 * t17 * t26 - 0.2e1 * FA0 * t20 * t26 + FH0 * t11 * t2 + FH0 * t17 * t2 + FH0 * t2 * t20 - 0.2e1 * FH0 * t2 * t3 - 0.2e1 * FH0 * t2 * t7;
  t43 = DenS(p1+p2, mH, GammaH);
  t44 = RE(t43);
  return(-0.32e2 * t39 / t1 / t3 * t44 * PREF_B_PHIxQCD);
}
