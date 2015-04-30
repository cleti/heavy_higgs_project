
#include "AMP_HEADER.h"

double Eval_B_QCDxQCD_GG (AMP_ARGS)
{

  
  AMP_DEFINITIONS;
  AP_REFS_B(ap);
  
  double t1;
  double t11;
  double t13;
  double t14;
  double t15;
  double t16;
  double t18;
  double t2;
  double t21;
  double t22;
  double t23;
  double t3;
  double t31;
  double t5;
  double t51;
  double t56;
  double t59;
  double t6;
  double t60;
  double t63;
  double t64;
  double t65;
  double t7;
  double t8;
  t1 = beta * beta;
  t2 = y * y;
  t3 = t1 * t2;
  t5 = sp(s1, s2);
  t6 = t1 * t1;
  t7 = t6 * t5;
  t8 = t2 * t2;
  t11 = t6 * s;
  t13 = sp(s1, p1);
  t14 = sp(s2, p2);
  t15 = t13 * t14;
  t16 = t1 * beta;
  t18 = t16 * t2 * y;
  t21 = sp(s1, p2);
  t22 = sp(s2, p1);
  t23 = t21 * t22;
  t31 = t16 * y;
  t51 = t1 * s;
  t56 = -0.2e1 * s * t2 * t7 + s * t7 * t8 - 0.4e1 * t1 * t13 * t14 - 0.4e1 * t1 * t21 * t22 - 0.2e1 * s * t3 + s * t5 + 0.2e1 * s * t7 + 0.2e1 * t11 * t2 - t11 * t8 + 0.4e1 * t15 * t18 + 0.4e1 * t15 * t3 - 0.4e1 * t15 * t31 - 0.4e1 * t18 * t23 + 0.4e1 * t3 * t23 + 0.4e1 * t23 * t31 - 0.2e1 * t5 * t51 + s - 0.2e1 * t11 + 0.2e1 * t51;
  t59 = pow(beta_y - 0.1e1, 0.2e1);
  t60 = 0.1e1 / t59;
  t63 = pow(beta_y + 0.1e1, 0.2e1);
  t64 = 0.1e1 / t63;
  t65 = 0.1e1 / s;
  return(0.16e2 * (t3 + 0.1e1) * t56 * t60 * t64 * t65 * PREF_B_QCDxQCD_CA - 0.8e1 * t56 * PREF_B_QCDxQCD_CF * t60 * t64 * t65);
}
