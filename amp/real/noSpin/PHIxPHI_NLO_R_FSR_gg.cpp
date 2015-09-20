
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_FSR (AMP_ARGS)
{
  AMP_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_R(ap);
  
  double t10;
  double t12;
  double t15;
  double t16;
  double t17;
  double t18;
  double t19;
  double t20;
  double t24;
  double t25;
  double t30;
  double t39;
  double t40;
  double t45;
  double t5;
  double t57;
  double t6;
  double t7;
  double t8;
  double t9;
  t5 = 0.1024e4 * At2_fA2_De + 0.256e3 * At2_fH2_De + 0.1024e4 * Bt2_fA2_De + 0.256e3 * Bt2_fH2_De;
  t6 = t5 * PREF_R_PHI_CF;
  t7 = sp(p1, p2);
  t8 = t7 * t7;
  t9 = sp(k2, p3);
  t10 = 0.1e1 / t9;
  t12 = sp(k1, p3);
  t15 = 0.2048e4 * At2_fA2_De;
  t16 = 0.512e3 * At2_fH2_De;
  t17 = 0.2048e4 * Bt2_fA2_De;
  t18 = 0.512e3 * Bt2_fH2_De;
  t19 = t15 + t16 + t17 + t18;
  t20 = t19 * PREF_R_PHI_CF;
  t24 = t8 * t7;
  t25 = -t19 * PREF_R_PHI_CF * t24;
  t30 = (0.4096e4 * At2_fA2_De + 0.1024e4 * At2_fH2_De) * PREF_R_PHI_CF * t8;
  t39 = -t5 * PREF_R_PHI_CF * t24 + (t15 + t16) * PREF_R_PHI_CF * t8;
  t40 = t9 * t9;
  t45 = t8 * t8;
  t57 = t12 * t12;
  return(t6 * t8 * t10 * t12 + t20 * t8 + (t25 + t30) * t10 + t39 / t40 + (t6 * t8 * t9 + t25 + t30 + (t20 * t45 + (-0.6144e4 * At2_fA2_De - 0.1536e4 * At2_fH2_De - t17 - t18) * PREF_R_PHI_CF * t24 + t30) * t10) / t12 + t39 / t57);
}
