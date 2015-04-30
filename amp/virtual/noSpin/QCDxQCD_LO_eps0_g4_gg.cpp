
#include "AMP_HEADER.h"

#include <math.h>

double Eval_B_QCDxQCD_GG (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_B(ap);


  double t1;
  double t12;
  double t16;
  double t19;
  double t2;
  double t21;
  double t3;
  double t4;
  double t9;
  t1 = beta * beta;
  t2 = t1 * t1;
  t3 = y * y;
  t4 = t3 * t3;
  t9 = t1 * t3;
  t12 = -0.2e1 * t2 * t3 + t2 * t4 - 0.2e1 * t1 + 0.2e1 * t2 + 0.2e1 * t9 - 0.1e1;
  t16 = pow(beta_y - 0.1e1, 0.2e1);
  t19 = pow(beta_y + 0.1e1, 0.2e1);
  t21 = 0.1e1 / t16 / t19;
  return (0.32e2 * t12 / Constants::CA * t21 * PREF_B_QCDxQCD - 0.16e2 * (t9 + 0.1e1) * t12 * t21 * PREF_B_QCDxQCD_CA);
}
