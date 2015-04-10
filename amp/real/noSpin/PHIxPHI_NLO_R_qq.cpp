
#include "AMP_HEADER.h"

double Eval_R_PHIxPHI_QQ (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t17;
  double t2;
  double t3;
  double t4;
  double t7;
  t1 = sp(p1, p3);
  t2 = t1 * t1;
  t3 = sp(p2, p3);
  t4 = t3 * t3;
  t7 = sp(k1, k2);
  t17 = sp(p1, p2);
  return(0.128e3 * PREF_V_CA * (t4 + t2) * (0.4e1 * t7 * At2_fA2_De + t7 * At2_fH2_De + 0.4e1 * t7 * Bt2_fA2_De + t7 * Bt2_fH2_De - 0.4e1 * At2_fA2_De - At2_fH2_De + 0.4e1 * Bt2_fA2_De + Bt2_fH2_De) / t17);
}
