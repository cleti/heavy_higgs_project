
#include "../../../../inc/Functions_Shared.h"

double Eval_B_2PHIxQCD (
  PS_2_2 const& ps)
{

#include "EVAL_V_PS_REFS"
  
  double t1;
  double t10;
  double t11;
  c_double t12;
  double t13;
  double t14;
  double t15;
  double t19;
  double t2;
  double t20;
  double t21;
  double t28;
  double t3;
  double t39;
  double t4;
  double t40;
  double t44;
  double t6;
  double t9;
  t1 = sp(k2, p1);
  t2 = sp(p1, p2);
  t3 = t2 * t1;
  t4 = sp(k1, p1);
  t6 = t4 * t4;
  t9 = 0.1e1 / t1;
  t10 = Bt * (t1 * t4 + t3 + t6) * t9;
  t11 = 0.1e1 / t4;
  t12 = DenS(p1+p2, mH, GammaH);
  t13 = RE(t12);
  t14 = t11 * t13;
  t15 = RE(F_ggH_p);
  t19 = IM(t12);
  t20 = t11 * t19;
  t21 = IM(F_ggH_p);
  t28 = t2 * t2;
  t39 = At * (t1 * t28 - t6 * t1 - 0.2e1 * t2 * t4 + 0.2e1 * t2 * t6 + t3 * t4 - t4 * t6 - 0.2e1 * t3) / t2 * t9;
  t40 = RE(F_ggH_s);
  t44 = IM(F_ggH_s);
  double ret = PREF_B_PHIxQCD*(0.4e1 * t10 * t14 * t15  - 0.4e1 * t10 * t20 * t21 + 0.4e1 * t14 * t39 * t40 - 0.4e1 * t20 * t39 * t44);
  return ret;
}
