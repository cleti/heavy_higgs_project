
#include "UID_HEADER.h"

double Eval_UID_PHIxPHI_SE00 (UID_ARGS)
{

UID_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_R(ap);


  double t11;
  double t12;
  double t19;
  double t3;
  double t4;
  double t6;
  double t7;
  double t9;
  t3 = sp(K1, K2);
  t4 = sp(P2, p1);
  t6 = sp(P2, p3);
  t7 = sp(p1, p3);
  t9 = t3 * t4 * t6 * t7;
  t11 = t4 * t6 * t7;
  t12 = t9 - t11;
  t19 = t9 + t11;
  return(0.4096e4 * t12 * At2_fA2_De * PREF_R_PHI_CA + 0.1024e4 * t12 * At2_fH2_De * PREF_R_PHI_CA + 0.4096e4 * t19 * Bt2_fA2_De * PREF_R_PHI_CA + 0.1024e4 * t19 * Bt2_fH2_De * PREF_R_PHI_CA);
}
