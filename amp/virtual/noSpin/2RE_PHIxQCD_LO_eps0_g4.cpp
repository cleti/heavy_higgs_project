
#include "AMP_HEADER.h"

double Eval_B_2PHIxQCD (AMP_ARGS)
{

  AMP_DEFINITIONS;
  HP_REFS_PHIxQCD(hp);
  AP_REFS_B(ap);

  double t7;
  t7 = 0.1e1 / (beta_y - 0.1e1) / (beta_y + 0.1e1) * PREF_B_PHIxQCD;
  return (-0.128e3 * Bt_fA_re * s * t7 - 0.64e2 * At_fH_re * (0.4e1 - s) * t7);
}
