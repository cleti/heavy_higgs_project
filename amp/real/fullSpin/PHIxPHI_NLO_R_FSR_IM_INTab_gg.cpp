
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_FSR_IM_INTab (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_R(ap);

  // need the imaginary parts of these prefactors here !!!
  double const& At_Bt_fH2_De = hp.At_Bt_fH2_DeIM;
  double const& At_Bt_fA2_De = hp.At_Bt_fA2_DeIM;

  double t11;
  double t12;
  double t16;
  double t18;
  double t19;
  double t23;
  double t24;
  double t25;
  double t29;
  double t3;
  double t33;
  double t39;
  double t4;
  double t40;
  double t43;
  double t5;
  double t6;
  double t63;
  double t7;
  double t71;
  double t8;
  t3 = 0.2048e4 * At_Bt_fA2_De + 0.512e3 * At_Bt_fH2_De;
  t4 = t3 * PREF_R_PHI_CF;
  t5 = sp(p1, p2);
  t6 = t5 * t5;
  t7 = sp(k2, p3);
  t8 = 0.1e1 / t7;
  t11 = t7 * t7;
  t12 = 0.1e1 / t11;
  t16 = sp(k1, s2);
  t18 = sp(s1, k2);
  t19 = t18 * t6;
  t23 = -t3 * PREF_R_PHI_CF;
  t24 = sp(s2, p3);
  t25 = t6 * t5;
  t29 = sp(s1, p3);
  t33 = (t18 * t4 + t24 * t4 + t29 * t4) * t6;
  t39 = -0.4096e4 * At_Bt_fA2_De - 0.1024e4 * At_Bt_fH2_De;
  t40 = t39 * PREF_R_PHI_CF;
  t43 = -t39 * PREF_R_PHI_CF;
  t63 = sp(k1, p3);
  t71 = t63 * t63;
  return((t12 * t4 * t6 + t4 * t6 * t8) * t16 + t4 * t19 * t8 + (t23 * t24 * t25 + t33) * t12 + ((t4 * t6 + (t25 * t40 + t43 * t6) * t8) * t16 + t4 * t19 + ((t18 * t40 + t24 * t23 + t23 * t29) * t25 + (t18 * t43 + t24 * t43 + t29 * t43) * t6) * t8) / t63 + (t16 * t4 * t6 + t23 * t25 * t29 + t33) / t71);
}
