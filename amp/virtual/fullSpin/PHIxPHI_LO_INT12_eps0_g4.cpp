
#include "AMP_HEADER.h"

double Eval_B_PHIxPHI_INT12 (AMP_ARGS)
{

AMP_DEFINITIONS
AMP_H12_REFS

  double t1;
  double t10;
  double t12;
  double t16;
  double t17;
  double t18;
  double t19;
  double t20;
  double t21;
  double t22;
  double t23;
  double t25;
  double t29;
  double t30;
  double t32;
  double t39;
  double t40;
  c_double t41;
  double t42;
  double t44;
  double t50;
  double t51;
  c_double t54;
  double t55;
  double t58;
  double t6;
  double t7;
  double t8;
  t1 = s * s;
  t6 = t1 * (0.4e1 * fA1 * fA2 + fH1 * fH2);
  t7 = At1 * At2;
  t8 = Bt1 * Bt2;
  t10 = sp(s1, k2);
  t12 = sp(k1, s2);
  t16 = beta * beta;
  t17 = t16 * s;
  t18 = t7 * t17;
  t19 = t8 * t17;
  t20 = t7 * s;
  t21 = t8 * s;
  t22 = 0.4e1 * t7;
  t23 = 0.4e1 * t8;
  t25 = sp(s1, s2);
  t29 = At1 * Bt2;
  t30 = At2 * Bt1;
  t32 = EPS_(k1, k2, s1, s2);
  t39 = 0.32e2 * t6 * (t7 - t8) * t10 * t12 - 0.8e1 * t6 * (t18 - t19 + t20 - t21 - t22 - t23) * t25 + 0.32e2 * t6 * (t29 + t30) * t32 + 0.8e1 * t6 * (t18 + t19 + t20 + t21 - t22 + t23);
  t40 = s;
  t41 = DenS(t40, mH2, GammaH2);
  t42 = RE(t41);
  t44 = t29 - t30;
  t50 = -0.32e2 * t10 * t44 * t6 - 0.32e2 * t12 * t44 * t6;
  t51 = IM(t41);
  t54 = DenS(t40, mH1, GammaH1);
  t55 = RE(t54);
  t58 = IM(t54);
  return(((t39 * t42 + t50 * t51) * t55 - t50 * t58 * t42 + t39 * t51 * t58) * PREF_B_PHIxPHI);
}
