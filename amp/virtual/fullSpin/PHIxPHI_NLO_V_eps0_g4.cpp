
#include "AMP_HEADER.h"

#include <math.h>

double Eval_V_PHIxPHI (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t10;
  double t11;
  double t123;
  double t130;
  double t14;
  double t18;
  double t19;
  double t20;
  double t24;
  double t25;
  double t26;
  double t27;
  double t3;
  double t41;
  double t42;
  double t43;
  double t5;
  double t51;
  double t6;
  double t64;
  double t70;
  double t73;
  double t76;
  double t8;
  double t87;
  double t9;
  double t91;
  double t99;
  t1 = s * s;
  t3 = s * Bt2_fH2_De;
  t5 = EPS_(k1, k2, s1, s2);
  t6 = At_Bt_fH2_De * t5;
  t8 = CA * Pi2;
  t9 = beta * beta;
  t10 = t9 * s;
  t11 = t10 * At2_fA2_De;
  t14 = t10 * At2_fH2_De;
  t18 = log(MUR2 / s);
  t19 = t18 * t18;
  t20 = CA * t19;
  t24 = sp(k1, s2);
  t25 = CA * t24;
  t26 = sp(s1, k2);
  t27 = t26 * t19;
  t41 = sp(s1, s2);
  t42 = CA * t41;
  t43 = t19 * s;
  t51 = t26 * Pi2;
  t64 = Pi2 * s;
  t70 = At_Bt_fA2_De * t5;
  t73 = t24 * t26;
  t76 = -0.8e1 * t25 * t51 * At2_fA2_De - 0.2e1 * t25 * t51 * At2_fH2_De + 0.8e1 * t25 * t51 * Bt2_fA2_De + 0.2e1 * t25 * t51 * Bt2_fH2_De + 0.4e1 * t42 * t43 * Bt2_fA2_De + t42 * t43 * Bt2_fH2_De - 0.4e1 * t42 * t64 * Bt2_fA2_De - t42 * t64 * Bt2_fH2_De - 0.64e2 * CA * t70 + 0.16e2 * t11 * t42 - 0.22e2 * t73 * At2_fH2_De;
  t87 = t42 * t19;
  t91 = t42 * Pi2;
  t99 = s * Bt2_fA2_De;
  t123 = CA * t5;
  t130 = 0.11e2 * s * t41 * t9 * At2_fH2_De - 0.16e2 * Pi2 * t123 * At_Bt_fA2_De - 0.4e1 * Pi2 * t123 * At_Bt_fH2_De - 0.32e2 * t25 * t26 * At2_fA2_De + 0.32e2 * t25 * t26 * Bt2_fA2_De + t20 * t3 + 0.4e1 * t20 * t6 + 0.16e2 * t20 * t70 + 0.4e1 * t20 * t99 - t3 * t8 - 0.16e2 * t42 * t99;
  return(-0.4e1 * t1 * PREF_V_PHI_CA * (0.8e1 * t25 * t27 * At2_fA2_De + 0.2e1 * t25 * t27 * At2_fH2_De - 0.8e1 * t25 * t27 * Bt2_fA2_De - 0.2e1 * t25 * t27 * Bt2_fH2_De + 0.4e1 * t20 * t11 - 0.4e1 * t8 * t11 + t20 * t14 - t8 * t14 - 0.11e2 * t3 - 0.44e2 * t6 + t76 - 0.16e2 * CA * t9 * s * At2_fA2_De - 0.16e2 * CA * s * Bt2_fA2_De - 0.11e2 * t41 * s * Bt2_fH2_De - 0.4e1 * t87 * t11 + 0.4e1 * t91 * t11 - t87 * t14 + t91 * t14 + 0.22e2 * t73 * Bt2_fH2_De - 0.4e1 * t8 * t99 + t130 - 0.11e2 * t14) / CA);
}
