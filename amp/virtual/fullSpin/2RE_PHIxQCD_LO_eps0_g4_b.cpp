
#include "AMP_HEADER.h"

double Eval_B_2PHIxQCD (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_B(ap);
  HP_REFS_PHIxQCD(hp);

  double t1;
  double t11;
  double t14;
  double t15;
  double t18;
  double t2;
  double t3;
  double t31;
  double t35;
  double t5;
  double t6;
  double t8;
  t1 = beta * beta;
  t2 = s * t1;
  t3 = sp(s1, s2);
  t5 = sp(s1, k2);
  t6 = sp(k1, s2);
  t8 = 0.2e1 * t5 * t6;
  t11 = y * y;
  t14 = 0.1e1 / (t1 * t11 - 0.1e1);
  t15 = PREF_B_PHIxQCD * t14;
  t18 = EPS_(k1, k2, s1, s2);
  t31 = t6 + t5;
  t35 = PREF_B_PHIxQCD * (t1 - 0.1e1) * t14;
  return(0.64e2 * At_fH_re * (-t3 * t2 + t2 + t8) * t15 - 0.256e3 * At_fA_re * t18 * t15 + 0.128e3 * Bt_fH_re * t18 * t15 + 0.128e3 * Bt_fA_re * (-s * t3 - s + t8) * t15 + 0.64e2 * At_fA_im * s * t31 * t35 + 0.32e2 * Bt_fH_im * s * t31 * t35);
}
