
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_FSR_IM_INTab (AMP_ARGS)
{

AMP_DEFINITIONS

  // need the imaginary parts of these prefactors here !!!
  double const& At_Bt_fH2_De = At_Bt_fH2_DeIM;
  double const& At_Bt_fA2_De = At_Bt_fA2_DeIM;

  double t1;
  double t11;
  double t13;
  double t14;
  double t15;
  double t19;
  double t2;
  double t21;
  double t23;
  double t29;
  double t3;
  double t37;
  double t4;
  double t49;
  double t5;
  double t7;
  double t8;
  t1 = sp(p1, p2);
  t2 = t1 * t1;
  t3 = sp(s1, k2);
  t4 = sp(k1, p3);
  t5 = t4 * t4;
  t7 = t3 * t4;
  t8 = sp(k2, p3);
  t11 = t8 * t8;
  t13 = sp(k1, s2);
  t14 = t13 * t5;
  t15 = t13 * t4;
  t19 = sp(s1, p3);
  t21 = sp(s2, p3);
  t23 = t4 * t8;
  t29 = -t11 * t13 - t11 * t19 - t11 * t3 - 0.2e1 * t15 * t8 - 0.2e1 * t19 * t23 - t19 * t5 - 0.2e1 * t21 * t23 - t21 * t5 - t3 * t5 - 0.2e1 * t7 * t8 - t14;
  t37 = t1 * t4;
  t49 = t1 * t11 * t19 + t1 * t13 * t23 + t1 * t21 * t5 + t1 * t23 * t3 + t19 * t37 * t8 + t21 * t37 * t8 - t3 * t5 * t8 - t11 * t15 - t11 * t21 - t11 * t7 - t14 * t8;
  return(-0.512e3 * t2 * (t29 + t49) * PREF_R_PHI_CF * (At_Bt_fH2_De + 0.4e1 * At_Bt_fA2_De) / t11 / t5);
}
