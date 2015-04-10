
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_QG (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t16;
  double t17;
  double t18;
  double t2;
  double t3;
  double t30;
  double t38;
  double t4;
  double t45;
  double t47;
  double t7;
  double t8;
  double t9;
  t1 = sp(p1, p2);
  t2 = t1 * t1;
  t3 = sp(p2, p3);
  t4 = t3 * t3;
  t7 = sp(k1, k2);
  t8 = sp(s1, s2);
  t9 = t7 * t8;
  t16 = sp(k1, s2);
  t17 = sp(s1, k2);
  t18 = t16 * t17;
  t30 = -0.4e1 * t18 * At2_fA2_De - t18 * At2_fH2_De + 0.4e1 * t18 * Bt2_fA2_De + t18 * Bt2_fH2_De - 0.4e1 * t7 * At2_fA2_De - t7 * At2_fH2_De - 0.4e1 * t7 * Bt2_fA2_De + 0.4e1 * t9 * At2_fA2_De + t9 * At2_fH2_De - 0.4e1 * t9 * Bt2_fA2_De - t9 * Bt2_fH2_De;
  t38 = EPS_(k1, k2, s1, s2);
  t45 = -0.8e1 * t38 * At_Bt_fA2_De - 0.2e1 * t38 * At_Bt_fH2_De - t7 * Bt2_fH2_De - 0.4e1 * t8 * At2_fA2_De - t8 * At2_fH2_De - 0.4e1 * t8 * Bt2_fA2_De - t8 * Bt2_fH2_De + 0.4e1 * At2_fA2_De + At2_fH2_De - 0.4e1 * Bt2_fA2_De - Bt2_fH2_De;
  t47 = sp(p1, p3);
  return(-0.128e3 * PREF_V_CA * (t2 + t4) * (t30 + t45) / t47);
}
