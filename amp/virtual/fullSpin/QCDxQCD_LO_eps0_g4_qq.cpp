
#include "AMP_HEADER.h"

double Eval_B_QCDxQCD_QQ (AMP_ARGS)
{

  AMP_DEFINITIONS;
  AP_REFS_B(ap);

  
  double t1;
  double t11;
  double t14;
  double t16;
  double t19;
  double t2;
  double t4;
  double t9;
  t1 = sp(s1, s2);
  t2 = beta * beta;
  t4 = y * y;
  t9 = sp(s1, p1);
  t11 = sp(s2, p2);
  t14 = sp(s2, p1);
  t16 = sp(s1, p2);
  t19 = t2 * s;
  return(-0.8e1 * (s * t1 * t2 * t4 + 0.4e1 * beta_y * t11 * t9 - 0.4e1 * beta_y * t14 * t16 - s * t2 * t4 - t1 * t19 + 0.4e1 * t11 * t9 + 0.4e1 * t14 * t16 - 0.2e1 * s + t19) * PREF_B_QCDxQCD / s);
}
