
#include "UID_HEADER.h"

double Eval_UID_PHIxPHI_Q_QG_ES00 (UID_ARGS)
{

  UID_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_R(ap);


  double t12;
  double t3;
  double t5;
  double t6;
  double t8;
  double t9;
  t3 = sp(p2, p3);
  t5 = sp(K1, p3);
  t6 = sp(K2, p3);
  t8 = sp(K1, p2);
  t9 = sp(K2, p2);
  t12 = sp(K1, K2);
  return (0.256e3 * t3 * PREF_R_PHI_CA / CA * (t5 - t3 + t6) * (t9 + t8) * (0.4e1 * t12 * At2_fA2_De + t12 * At2_fH2_De + 0.4e1 * t12 * Bt2_fA2_De + t12 * Bt2_fH2_De - 0.4e1 * At2_fA2_De - At2_fH2_De + 0.4e1 * Bt2_fA2_De + Bt2_fH2_De));
}
