
#include "UID_HEADER.h"

double Eval_UID_PHIxPHI_Q_QG_SE00 (UID_ARGS)
{

  UID_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_R(ap);


  double t15;
  double t16;
  double t19;
  double t20;
  double t3;
  double t5;
  t3 = sp(p1, p3);
  t5 = sp(K1, K2);
  t15 = sp(K1, p3);
  t16 = sp(K2, p3);
  t19 = sp(K1, p1);
  t20 = sp(K2, p1);
  return(0.256e3 * t3 * PREF_R_PHI_CA / CA * (0.4e1 * t5 * At2_fA2_De + t5 * At2_fH2_De + 0.4e1 * t5 * Bt2_fA2_De + t5 * Bt2_fH2_De - 0.4e1 * At2_fA2_De - At2_fH2_De + 0.4e1 * Bt2_fA2_De + Bt2_fH2_De) * (t15 + t16 - t3) * (t19 + t20));
}
