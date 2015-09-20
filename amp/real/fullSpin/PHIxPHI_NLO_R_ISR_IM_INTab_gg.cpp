
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_ISR_IM_INTab (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_R(ap);

  // need the imaginary parts of these prefactors here !!!
  double const& At_Bt_fH2_De = hp.At_Bt_fH2_DeIM;
  double const& At_Bt_fA2_De = hp.At_Bt_fA2_DeIM;

  double t1;
  double t11;
  double t17;
  double t2;
  double t20;
  double t29;
  double t3;
  double t32;
  double t39;
  double t4;
  double t40;
  double t41;
  double t42;
  double t5;
  double t8;
  t1 = sp(p1, p3);
  t2 = t1 * t1;
  t3 = t2 * t2;
  t4 = sp(p2, p3);
  t5 = t2 * t1;
  t8 = sp(p1, p2);
  t11 = t4 * t4;
  t17 = t8 * t8;
  t20 = t11 * t4;
  t29 = t17 * t8;
  t32 = t11 * t11;
  t39 = t17 * t17;
  t40 = -0.6e1 * t1 * t11 * t8 + 0.6e1 * t1 * t17 * t4 - 0.6e1 * t2 * t4 * t8 + 0.2e1 * t1 * t20 - 0.2e1 * t1 * t29 + 0.3e1 * t11 * t17 + 0.3e1 * t11 * t2 + 0.3e1 * t17 * t2 - 0.2e1 * t20 * t8 - 0.2e1 * t29 * t4 + 0.2e1 * t4 * t5 - 0.2e1 * t5 * t8 + t3 + t32 + t39;
  t41 = sp(s1, k2);
  t42 = sp(k1, s2);
  return(-0.1024e4 * t40 * (t41 + t42) * PREF_R_PHI_CA * (At_Bt_fH2_De + 0.4e1 * At_Bt_fA2_De) / t4 / t1 / t8);
}
