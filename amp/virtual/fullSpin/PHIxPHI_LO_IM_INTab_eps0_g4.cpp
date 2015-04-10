
// this is the imaginary part of the At-Bt interference of a single Higgs boson

#include "AMP_HEADER.h"

double Eval_B_PHIxPHI_IM_INTab (AMP_ARGS)
{

  AMP_DEFINITIONS
  
  double t2;
  double t3;
  double t5;
  double t7;
  t2 = s * s;
  t3 = sp(k1, s2);
  t5 = sp(s1, k2);
  t7 = -t2 * (t3 + t5);
  return 32.0 * (4.0 * At_Bt_fA2_DeIM + At_Bt_fH2_DeIM ) * t7 * PREF_B_PHIxPHI;
}
