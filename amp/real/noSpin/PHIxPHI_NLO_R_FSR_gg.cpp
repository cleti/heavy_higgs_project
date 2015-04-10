
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_FSR (AMP_ARGS)
{

AMP_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_R(ap);


  double t1;
  double t10;
  double t102;
  double t109;
  double t11;
  double t114;
  double t12;
  double t14;
  double t16;
  double t17;
  double t2;
  double t20;
  double t23;
  double t26;
  double t29;
  double t36;
  double t39;
  double t4;
  double t5;
  double t59;
  double t65;
  double t68;
  double t78;
  double t88;
  double t92;
  double t97;
  t1 = sp(p1, p2);
  t2 = t1 * t1;
  t4 = sp(k2, p3);
  t5 = t4 * t4;
  t10 = sp(k1, p3);
  t11 = t10 * t10;
  t12 = t11 * At2_fA2_De;
  t14 = t11 * At2_fH2_De;
  t16 = t5 * t1;
  t17 = t10 * At2_fA2_De;
  t20 = t10 * At2_fH2_De;
  t23 = t10 * Bt2_fA2_De;
  t26 = t10 * Bt2_fH2_De;
  t29 = t4 * t2;
  t36 = -0.8e1 * t16 * t17 - 0.2e1 * t16 * t20 - 0.8e1 * t16 * t23 - 0.2e1 * t16 * t26 + 0.8e1 * t17 * t29 + 0.2e1 * t20 * t29 + 0.8e1 * t23 * t29 + 0.8e1 * t5 * At2_fA2_De + 0.2e1 * t5 * At2_fH2_De + 0.8e1 * t12 + 0.2e1 * t14;
  t39 = t4 * t1;
  t59 = t5 * t4 * t10;
  t65 = -0.8e1 * t11 * t39 * Bt2_fA2_De - 0.2e1 * t11 * t39 * Bt2_fH2_De - 0.8e1 * t12 * t39 - 0.2e1 * t14 * t39 - 0.24e2 * t17 * t39 - 0.6e1 * t20 * t39 - 0.8e1 * t23 * t39 + 0.2e1 * t26 * t29 - 0.2e1 * t26 * t39 + 0.4e1 * t59 * At2_fA2_De + t59 * At2_fH2_De + 0.4e1 * t59 * Bt2_fA2_De;
  t68 = t5 * t11;
  t78 = t4 * t11 * t10;
  t88 = -0.4e1 * t16 * At2_fA2_De - t16 * At2_fH2_De + t59 * Bt2_fH2_De + 0.8e1 * t68 * At2_fA2_De + 0.2e1 * t68 * At2_fH2_De + 0.8e1 * t68 * Bt2_fA2_De + 0.2e1 * t68 * Bt2_fH2_De + 0.4e1 * t78 * At2_fA2_De + t78 * At2_fH2_De + 0.4e1 * t78 * Bt2_fA2_De + t78 * Bt2_fH2_De;
  t92 = t5 * t10;
  t97 = t4 * t11;
  t102 = t1 * t11;
  t109 = t4 * t10;
  t114 = -0.4e1 * t102 * At2_fA2_De - t102 * At2_fH2_De - 0.4e1 * t102 * Bt2_fA2_De - t102 * Bt2_fH2_De + 0.16e2 * t109 * At2_fA2_De + 0.4e1 * t109 * At2_fH2_De - 0.4e1 * t16 * Bt2_fA2_De - t16 * Bt2_fH2_De + 0.16e2 * t92 * At2_fA2_De + 0.4e1 * t92 * At2_fH2_De + 0.16e2 * t97 * At2_fA2_De + 0.4e1 * t97 * At2_fH2_De;
  return(0.256e3 * t2 * PREF_R_PHI_CF * (t36 + t65 + t88 + t114) / t5 / t11);
}
