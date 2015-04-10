
#include "UID_HEADER.h"

double Eval_UID_PHIxPHI_IM_INTab_SE00 (UID_ARGS)
{

UID_DEFINITIONS


  double t10;
  double t13;
  double t3;
  double t4;
  double t5;
  double t6;
  double t7;
  t3 = sp(p3, P2);
  t4 = sp(p1, p3);
  t5 = t3 * t4;
  t6 = sp(p1, P2);
  t7 = sp(K1, S2);
  t10 = sp(K2, S1);
  t13 = -t10 * t5 * t6 - t5 * t6 * t7;
  return(0.8192e4 * t13 * At_Bt_fA2_De * PREF_R_PHI_CA + 0.2048e4 * t13 * At_Bt_fH2_De * PREF_R_PHI_CA);
}
