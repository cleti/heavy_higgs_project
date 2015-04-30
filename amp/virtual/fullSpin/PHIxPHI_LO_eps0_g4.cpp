
#include "AMP_HEADER.h"

double Eval_B_PHIxPHI (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_B(ap);
  HP_REFS_PHIxPHI(hp);
  
  double t1;
  double t11;
  double t12;
  double t13;
  double t19;
  double t2;
  double t24;
  double t25;
  double t3;
  double t31;
  double t4;
  double t5;
  double t6;
  double t9;
  t1 = s * s;
  t2 = t1 * s;
  t3 = beta * beta;
  t4 = t2 * t3;
  t5 = sp(s1, s2);
  t6 = t4 * t5;
  t9 = sp(s1, k2);
  t11 = sp(k1, s2);
  t12 = t1 * t9 * t11;
  t13 = 0.16e2 * t12;
  t19 = 0.64e2 * t12;
  t24 = EPS_(k1, k2, s1, s2);
  t25 = t24 * PREF_B_PHIxPHI;
  t31 = t2 * t5;
  return(At2_fH2_De * (-0.8e1 * t6 + 0.8e1 * t4 + t13) * PREF_B_PHIxPHI + At2_fA2_De * (-0.32e2 * t6 + 0.32e2 * t4 + t19) * PREF_B_PHIxPHI + 0.32e2 * At_Bt_fH2_De * t1 * t25 + 0.128e3 * At_Bt_fA2_De * t1 * t25 + Bt2_fH2_De * (-t13 + 0.8e1 * t31 + 0.8e1 * t2) * PREF_B_PHIxPHI + Bt2_fA2_De * (-t19 + 0.32e2 * t31 + 0.32e2 * t2) * PREF_B_PHIxPHI);
}
