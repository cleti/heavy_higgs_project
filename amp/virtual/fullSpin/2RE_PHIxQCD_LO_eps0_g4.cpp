
#include "AMP_HEADER.h"

double Eval_B_2PHIxQCD (AMP_ARGS)
{

AMP_DEFINITIONS


  double t11;
  double t12;
  double t13;
  double t15;
  double t22;
  double t26;
  double t27;
  double t3;
  double t31;
  double t34;
  double t4;
  double t47;
  double t57;
  double t6;
  double t7;
  double t8;
  t3 = 0.1e1 / (beta_y - 0.1e1);
  t4 = (0.4e1 - s) * t3;
  t6 = 0.1e1 / (beta_y + 0.1e1);
  t7 = sp(s1, s2);
  t8 = t6 * t7;
  t11 = t3 * t6;
  t12 = sp(s1, k2);
  t13 = sp(k1, s2);
  t15 = t11 * t12 * t13;
  t22 = beta * beta;
  t26 = (s * t22 - s + 0.4e1) / s;
  t27 = EPS_(k2, p1, s1, s2);
  t31 = EPS_(k1, k2, s1, s2);
  t34 = EPS_(k1, p1, s1, s2);
  t47 = s * t3;
  t57 = -t11 * t12 - t11 * t13;
  return(At_fH_re * (-0.64e2 * t4 * t6 + 0.64e2 * t4 * t8 + 0.128e3 * t15) * PREF_B_PHIxQCD + At_fA_re * (-0.128e3 * t11 * t26 * t27 + 0.128e3 * t11 * t26 * t34 - 0.256e3 * t11 * t31) * PREF_B_PHIxQCD + 0.128e3 * Bt_fH_re * t3 * t6 * t31 * PREF_B_PHIxQCD + Bt_fA_re * (-0.128e3 * t47 * t6 - 0.128e3 * t47 * t8 + 0.256e3 * t15) * PREF_B_PHIxQCD + 0.256e3 * At_fA_im * t57 * PREF_B_PHIxQCD + 0.128e3 * Bt_fH_im * t57 * PREF_B_PHIxQCD);
}
