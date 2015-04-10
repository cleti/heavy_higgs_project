
#include "AMP_HEADER.h"

#include <math.h>

double Eval_B_QCDxQCD_GG (AMP_ARGS)
{

  AMP_DEFINITIONS
  AP_REFS_B(ap)

  double t1;
  double t10;
  double t15;
  double t16;
  double t18;
  double t19;
  double t20;
  double t24;
  double t30;
  double t4;
  double t5;
  double t6;
  double t7;
  double t8;
  double t9;
  t1 = beta * y;
  t4 = (t1 - 0.1e1) * (t1 + 0.1e1);
  t5 = beta * beta;
  t6 = t5 * t5;
  t7 = y * y;
  t8 = t7 * t7;
  t9 = t6 * t8;
  t10 = t6 * t7;
  t15 = pow(beta_y - 0.1e1, 0.2e1);
  t16 = 0.1e1 / t15;
  t18 = pow(beta_y + 0.1e1, 0.2e1);
  t19 = 0.1e1 / t18;
  t20 = t16 * t19;
  t24 = t6 * t5;
  t30 = t5 * t7;
  return(0.16e2 * t4 * (t9 - 0.2e1 * t10 + 0.1e1) * t20 * PREF_B_QCDxQCD_CA - 0.64e2 * (t24 * t7 * t8 + t24 * t7 - 0.2e1 * t24 * t8 - t10 + t30 - t5 + t6 + t9 - 0.1e1) * t16 * t19 * PREF_B_QCDxQCD_CF + 0.64e2 * t5 * (y - 0.1e1) * (y + 0.1e1) * (t30 - t5 + 0.1e1) * t4 * t20 * PREF_B_QCDxQCD_CFCA);
}
