
#include "AMP_HEADER.h"

double Eval_B_PHIxPHI (AMP_ARGS)
{

  AMP_DEFINITIONS;
  HP_REFS_PHIxPHI(hp);
  AP_REFS_B(ap);
    
  double t1;
  double t2;
  double t4;
  double t5;
  t1 = s * s;
  t2 = s * t1;
  t4 = beta * beta;
  t5 = t4 * PREF_B_PHIxPHI;
  return (0.32e2 * t2 * t5 * At2_fA2_De + 0.8e1 * t2 * t5 * At2_fH2_De + 0.32e2 * t2 * Bt2_fA2_De * PREF_B_PHIxPHI + 0.8e1 * t2 * Bt2_fH2_De * PREF_B_PHIxPHI);
}
