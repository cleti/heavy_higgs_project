
#include "AMP_HEADER.h"

double Eval_B_QCDxQCD_QQ (AMP_ARGS)
{

AMP_DEFINITIONS


  double t1;
  double t11;
  double t12;
  double t16;
  double t20;
  double t24;
  double t4;
  double t5;
  t1 = 0.1e1 / s;
  t4 = sp(s1, p1);
  t5 = sp(s2, p2);
  t11 = sp(s2, p1);
  t12 = sp(s1, p2);
  t16 = beta * beta;
  t20 = sp(s1, s2);
  t24 = y * y;
  return((-0.32e2 * t1 * (beta_y + 0.1e1) * t4 * t5 + 0.32e2 * (beta_y - 0.1e1) * t1 * t11 * t12 - 0.8e1 * t16 * (y - 0.1e1) * (y + 0.1e1) * t20 + 0.8e1 * t16 * t24 - 0.8e1 * t16 + 0.16e2) * PREF_B_QCDxQCD);
}
